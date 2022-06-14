//File that has the main function
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "spline.h"
#include "functions.h"
#include "particle.h"
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include "opacities.h"

using namespace std;
using std::fill;

vector <Particle> particles; //initiate vector of "Particle" (Object defined in particle.h)

//**********************************************************************
//**********************************************************************

//variable declaration for input file reading

int in_c = 0;
int tau_type, outflow;
int cont = 0;
string T_int_s;
string outflow_s;
string output_file;
double t_star, m_star, r_star;
double period, planet_mass, planet_radius;
double s_0;
double rmin, rmax, tmin, tmax, pmin,pmax;
double mdot_read, major_timestep, no_orbits, nparticles;
bool tau_constant;
string line;
ifstream input("dusty_tails.in", ios::out);
Opacities opac;

//**********************************************************************
//**********************************************************************

// variables declarations
//frame of reference parameters:
double ang_vel;
//Stellar parameters:
 double Mstar_kg, Mstar_sun, Rstar, Temp, lum;

//Planetary parameters:
 double Period_days, T, a, m_planet, r_planet, r_start, r_planet_dim, r_h, Temp_p;

//Dust parameters:
 double A, Bp, rho_d, mu, alpha;

//Outflow parameters:
 double mdot, v_esc, v_th, mu_gas;

//Some dimensionless quantitites:
 double G_dim, c_dim;

//Grid parameters:
double d_r_min, d_r_max,  d_t_min, d_t_max, d_p_min, d_p_max;

//variables for opacity:
string composition, comp;
string opac_data, opacity_dir;

//**********************************************************************
//**********************************************************************

//Variable declarations for optical depth grid building
//Cell numbers defined in constants.h

double r_a [r_cells + 1];
double r_b [r_cells];
double theta_a [t_cells + 1];
double theta_b [t_cells];
double phi_a [p_cells + 1];
double phi_b [p_cells];

double r_min = 0.0, r_max = r_cells_d;
double theta_min = 0.0, theta_max =  t_cells_d;
double phi_min = 0.0, phi_max = p_cells_d;

double dr = (r_max - r_min)/ r_cells_d;
double dtheta = (theta_max - theta_min ) / t_cells_d;
double dphi = ( phi_max - phi_min) / p_cells_d;

double*** extinction;
double*** optical_depth;

//**********************************************************************
//**********************************************************************

int main() {
  //Input file reading and needed variable declaration
  
  if (input.is_open())
  {
    while ( getline (input,line) )
    {
      if (in_c == 0) {
        cout << "Reading stellar parameters... " << endl;
        t_star = stod(line.substr(10,16));
        m_star = stod(line.substr(21,24));
        r_star = stod(line.substr(31,34));
        cout << "Tstar= " << t_star << "K, Mstar= " << m_star << " Msun, Rstar= "<< r_star << " Rsun." << endl;
        cout << "\n" ;
      }
      if (in_c == 1) {
        cout << "Reading planetary parameters... " << endl;
        period = stod(line.substr(11,15));
        planet_mass = stod(line.substr(22,25));
        planet_radius = stod(line.substr(35,38));
         cout << "Period= " << period << "hours, Mplanet= " << planet_mass << " Mearth, Rplanet= "<< planet_radius << " Rearth." << endl;
         cout << "\n" ;
      }
      if (in_c == 2) {
        cout << "Reading dust parameters... " << endl;
        s_0 = stod(line.substr(12,18)) * 1.0e-4;
        composition = line.substr(25,36);
        //cout << "Dust is composed of " << composition << endl;
        //cout << "initial dust grains size " << s_0 << " micron" << endl;
        
      }
      if (in_c == 3) {
        rmin = stod(line.substr(14,18));
        rmax = stod(line.substr(27,31));
        tmin = stod(line.substr(40,44));
        tmax = stod(line.substr(52,56));
        pmin = stod(line.substr(65,69)) * -1.0;
        pmax = stod(line.substr(78,82));

      }
      if (in_c == 4) {
        outflow = stoi(line.substr(14));
        mdot_read = stod(line.substr(23,26));
        tau_type = stoi(line.substr(36));
        if (outflow==1)  {
          cout << "Outflow is radially outwards from the whole planetary surface." << endl;
          outflow_s = "sph";
        } else {
          cout <<  "Outflow is from the planet's dayside." << endl;
          outflow_s = "dayside";
        }
        cout << "\n" ;
        cout << "The planetary mass loss rate is " << mdot_read << " Mearth/Gyr." << endl;
        cout << "\n" ;
        if (tau_type ==0) {
          cout << "Optical depth of dust is traced." << endl;
          tau_constant = false;
        } else {
          cout << "Optical depth of dust is constant and thin (0.1)." << endl;
          tau_constant = true;
        }
        cout << "\n" ;
      }
      if (in_c == 5) {
        major_timestep = stod(line.substr(10,14));
        cout << "DELTA T " << major_timestep << endl;
        no_orbits = stod(line.substr(25,28));
        nparticles = stoi(line.substr(38,41));
        cont = stoi(line.substr(50,51));
        cout << cont << endl;
      }

     
      
      in_c = in_c + 1;
      
    }
    input.close();
  }

//**********************************************************************
//**********************************************************************
// variables defined

//Stellar parameters:
 Mstar_kg = m_star*Msun; //mass of star in kg
 Mstar_sun = m_star; //mass of star in terms of mass of the sun
 Rstar = r_star; //stellar radius in sun radii
 Temp = t_star; //stars temperature
 cout << Temp << endl;
 lum = sigma*4.0*PI* pow(Rstar*Rsun_cgs, 2.0) * pow(Temp, 4.0); //stellar luminosity

//Planetary parameters:
 Period_days = period/24.0; //period of planet in days
 T = period*60.0*60.0; //period of planet in seconds
 a = pow((G*Mstar_kg* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in meters
 m_planet = (planet_mass*Mearth)/Mstar_kg; //planetary mass in solar masses
 r_planet = planet_radius*Rearth; //planetary radius in meters
 r_start = (2.*r_planet)/a; //start position for particles
 r_planet_dim = r_planet/a; //planetary radius in terms of semi major axis
 r_h = pow(m_planet/3.0, 1.0/3.0); //hill radius

// Frame of reference parameters:


//Outflow parameters:
 mdot =  mdot_read*Mearth_cgs/gyr;
 cout << "mdot " << mdot << endl;
 v_esc = (pow((2.0*G *planet_mass*Mearth)/(planet_radius*Rearth), 0.5)) * (T/a); //escape velocity
 if (outflow == 1) {
   //planet is not tidally locked
   Temp_p = pow(Rstar*Rsun_cgs / (2.0*a*100.0), 0.5) * Temp;  
 } else{
   //planet tidally locked
   Temp_p = pow((Rstar*Rsun_cgs) / (a*100.0), 0.5) * Temp;  
 }
 //gas is composed of a mixture of SiO, Mg, O, O2, Fe, SiO2 and MgO - with the fractions as indicated in Booth et al. 2022
 mu_gas = 0.281*60.083 + 0.250*24.305 + 0.223*15.999 + 0.158*32.0 + 0.079*55.845 + 0.005*60.08 + 0.003*40.3044;
 cout << "the planets temperature is " << Temp_p << endl;
 cout << "thermal vel" << sqrt((kb*Temp_p)/(mu_gas*amu)) << " cm/s " << endl;
 v_th = sqrt((kb*Temp_p)/(mu_gas*amu)) / (a*100);

cout << "thermal velocity in a_ps is " << v_th << " ap/s" << endl;
//Some dimensionless quantitites:
 G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant
 c_dim = clight * (T / a); //dimensionless speed of light
 ang_vel = pow((G_dim *(m_planet + 1.0)), 0.5);
//Grid parameters:
d_r_min = rmin;
d_r_max = rmax;
d_t_min = tmin;
d_t_max = tmax;
d_p_min = pmin;
d_p_max = pmax;

//Dust parameters
if (composition.substr(0,5) == "Al2O3") {
  cout << "Dust is composed of Corundum." << endl;
  opac_data = "corundum_K95";
  comp = "Al2O3";
  A = 7.74e+4; //clausius claperyon relation
  Bp = 39.3; //clausius claperyon relation
  rho_d = 4.0; //dust density
  mu = 101.961;
  alpha = 0.1;

} else if (composition.substr(0,7) == "Fe2SiO4") {
  cout << "Dust is composed of Fayalite." << endl;
  opac_data = "fayalite_F01";
  comp = "Fe2SiO4";
  A = 6.04e+4;
  Bp = 38.1;
  rho_d = 4.39;
  mu =  203.774;
  alpha = 0.1;

} else if (composition.substr(0,1)=="C") {
  cout << "Dust is composed of Graphite." << endl;
  opac_data = "graphite_D84";
  comp = "C";
  A = 9.36e+4;
  Bp = 36.2;
  rho_d = 2.16;
  mu =  12.011;
  alpha = 0.1;

} else if (composition.substr(0,6)=="MgSiO3"){
  cout << "Dust is composed of Enstatite." << endl;
  opac_data = "enstatite_J98_J94_D95";
  comp = "MgSiO3";
  A = 6.89e+4;
  Bp = 37.8;
  rho_d = 3.20;
  mu = 100.389;
  alpha = 0.1;

} else if (composition.substr(0,7)=="Mg2SiO4") {
  cout << "Dust is composed of Olivine." << endl;
  opac_data = "olivine_F01";
  comp = "Mg2SiO4";
  A = 6.53e+4;
  Bp = 34.3;
  rho_d = 3.27;
  mu = 140.694;
  alpha = 0.1;

} else if (composition.substr(0,3)=="SiC") {
  cout << "Dust is composed of silicon carbide." << endl;
  opac_data = "silicon_carbide_L93";
  comp = "SiC";
  A= 7.85e+4;
  Bp = 37.4;
  rho_d = 3.22;
  mu = 44.085;
  alpha = 0.04;


} else{
  cout << "Composition unknown, stopping.";
  abort();
}

cout << "rho= " << rho_d << " g/cm3, initial size= " << s_0 << " cm." << endl;
cout << "Clausius-Claperyon parameters, A= " << A << "K, B= " << Bp << endl;
cout << "\n" ;

opacity_dir = "./opacs_jankovic/calc_dust_opac/"+opac_data+"/opac_";


int T_int { static_cast<int> (Temp)};
T_int_s = to_string(T_int);

opac.read_data((opacity_dir+"temp.dat").c_str(), (opacity_dir+"sizes.dat").c_str(),
        (opacity_dir+"planck_abs_tstar"+T_int_s+".dat").c_str(), 
        (opacity_dir+"planck_sca_tstar"+T_int_s+".dat").c_str(),
        (opacity_dir+"planck_abs.dat").c_str(), (opacity_dir+"planck_sca.dat").c_str(),
        true);

        
//initiation parameters
long int total_particles = nparticles; //initial number of particles to start simulation with
double t_common = major_timestep;
double big_step = major_timestep; //big time step (in terms of planetary orbits)
double end_t = no_orbits; // end time of simulation
double total_t = 0.0; // total time that has passed, so 0 in the beginning
long int current_particles = 0; // number of current particles in simulation

//PROGRAM INITIATION
//grid built if optical depth is to be traced
  
  cout << "dr " << dr << endl;
  build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min); //grid for optical depth calculations
  //cout << "Built grid for optical depth tracing." << endl;
  

  //add first particles
  add_particles(particles, current_particles, total_particles, 0.0); // call function that adds particles to simulation (in particles.cpp file)
  cout << "cont " << cont << endl;
  cout << "current "<< current_particles << endl;
  cout << "total "<< total_particles << endl;
  //Solve particles is the main routine of the program. It is in particles.cpp.
  solve_particles(0.00, end_t, particles, total_particles,current_particles);
  //cout << planet_x << endl;
  
  //cout << opac.stellar_abs(0.2e-4) + opac.stellar_scat(0.2e-4) << endl;
  return 0;

}
