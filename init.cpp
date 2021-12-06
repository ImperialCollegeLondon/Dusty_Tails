//File that has the main function

#include <iostream>
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


using namespace std;
using std::fill;

vector <Particle> particles; //initiate vector of "Particle" (Object defined in particle.h)

//**********************************************************************
//**********************************************************************

//variable declaration for input file reading

int in_c = 0;
int tau_type, outflow;
double t_star, m_star, r_star;
double period, planet_mass, planet_radius;
double A_dust, B_dust, alpha_dust, rho_dust, s_0, mu_dust;
double rmin, rmax, tmin, tmax, pmin,pmax;
double mdot_read, major_timestep, no_orbits, nparticles;
bool tau_constant;
string line;
ifstream input("dusty_tails.in", ios::out);


//**********************************************************************
//**********************************************************************

// variables declarations
//Stellar parameters:
 double Mstar_kg, Mstar_sun, Rstar, Temp, lum;

//Planetary parameters:
 double Period_days, T, a, m_planet, r_planet, r_start, r_planet_dim, r_h;

//Dust parameters:
 double A, B, alpha, mu, rho_d;

//Outflow parameters:
 double mdot, v_esc;

//Some dimensionless quantitites:
 double G_dim, c_dim;

//Grid parameters:
double d_r_min, d_r_max,  d_t_min, d_t_max, d_p_min, d_p_max;

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

double extinction [r_cells][t_cells][p_cells];
double optical_depth [r_cells][t_cells][p_cells];

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
        A_dust = stod(line.substr(11,17));
        B_dust = stod(line.substr(22,26));
        alpha_dust = stod(line.substr(38,42));
        rho_dust = stod(line.substr(51,54));
        s_0 = stod(line.substr(65,70)) * 1.0e-4;
        mu_dust = stod(line.substr(79,86));
        cout << "rho= " << rho_dust << "g/cm3, initial size= " << s_0 << " cm." << endl;
        cout << "Sublimation coefficient= " << alpha_dust << ", molecular weight " << mu_dust << " gm" <<endl;
        cout << "Clausius-Claperyon parameters, A= " << A_dust << ", B= " << B_dust << endl;
        cout << "\n" ;
        
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
        } else {
          cout <<  "Outflow is from the planet's dayside." << endl;
        }
        cout << "\n" ;
        cout << "The planetary mass loss rate is " << mdot_read << " Mearth/Gyr." << endl;
        cout << "\n" ;
        if (tau_type ==1) {
          cout << "Optical depth of dust is traced." << endl;
          bool tau_constant = false;
        } else {
          cout << "Optical depth of dust is constant and thin (0.1)." << endl;
          bool tau_constant = true;
        }
        cout << "\n" ;
      }
      if (in_c == 5) {
        major_timestep = stod(line.substr(10,14));
        no_orbits = stod(line.substr(25,28));
        nparticles = stoi(line.substr(38,41));
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
 lum = sigma*4.0*PI* pow(Rstar*Rsun_cgs, 2.0) * pow(Temp, 4.0); //stellar luminosity

//Planetary parameters:
 Period_days = period/24.0; //period of planet in days
 T = period*60.0*60.0; //period of planet in seconds
 a = pow((G*Mstar_kg* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in meters
 m_planet = (planet_mass*Mearth)/Mstar_kg; //planetary mass in kg
 r_planet = planet_radius*Rearth; //planetary radius in meters
 r_start = (2.*r_planet)/a; //start position for particles
 r_planet_dim = r_planet/a; //planetary radius in terms of semi major axis
 r_h = pow(m_planet/3.0, 1.0/3.0); //hill radius

//Dust parameters:
 A = A_dust; //clausius claperyon relation
 B = B_dust; //clausius claperyon relation
 alpha = alpha_dust;
 mu = mu_dust; //molecular weight
 rho_d = rho_dust; //dust density

//Outflow parameters:
 mdot =  mdot_read*Mearth_cgs/gyr;
 v_esc = (pow((2.0*G *m_planet*Mearth)/(r_planet*Rearth), 0.5)) * (T/a); //escape velocity m/s

//Some dimensionless quantitites:
 G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant
 c_dim = clight * (T / a); //dimensionless speed of light

//Grid parameters:
d_r_min = rmin;
d_r_max = rmax;
d_t_min = tmin;
d_t_max = tmax;
d_p_min = pmin;
d_p_max = pmax;

//initiation parameters
long int total_particles = nparticles; //initial number of particles to start simulation with
double t_common = major_timestep;
double big_step = major_timestep; //big time step (in terms of planetary orbits)
double end_t = no_orbits; // end time of simulation
double total_t = 0.0; // total time that has passed, so 0 in the beginning
long int current_particles = 0; // number of current particles in simulation

//PROGRAM INITIATION
//grid built if optical depth is to be traced
  
  if (tau_constant == false) {
  build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min); //grid for optical depth calculations
  cout << "Built grid for optical depth tracing." << endl;
  }

  //add first particles
  add_particles(particles, current_particles, total_particles, 0.0); // call function that adds particles to simulation (in particles.cpp file)
  
  //Solve particles is the main routine of the program. It is in particles.cpp.
  solve_particles(0.00, end_t, particles, total_particles,current_particles);



  return 0;

}
