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

//Variable declarations for optical depth grid building
//Cell numbers defined in constants.h

//**********************************************************************
//**********************************************************************

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
  string line;
  ifstream input("dusty_tails.in", ios::out);
  int in_c = 0;
  int tau_type, outflow;
  double t_star, m_star, r_star;
  double period, planet_mass, planet_radius;
  double A_dust, B_dust, alpha_dust, rho_dust, s_0;
  double rmin, rmax, tmin, tmax, pmin,pmax;
  double mdot_read, major_timestep, no_orbits, nparticles;

  if (input.is_open())
  {
    while ( getline (input,line) )
    {
      if (in_c == 0) {
        cout << "Reading stellar parameters... " << endl;
        t_star = stod(line.substr(10,16));
        //cout<< t_star << endl;
        m_star = stod(line.substr(21,24));
        //cout<< m_star << endl;
        r_star = stod(line.substr(31,34));
        //cout<< r_star << endl;
        cout << "Tstar= " << t_star << "K, Mstar= " << m_star << " Msun, Rstar= "<< r_star << " Rsun." << endl;
      }
      if (in_c == 1) {
        cout << "Reading planetary parameters... " << endl;
        period = stod(line.substr(11,15));
        planet_mass = stod(line.substr(22,25));
        planet_radius = stod(line.substr(35,38));
         cout << "Period= " << period << "hours, Mplanet= " << planet_mass << " Mearth, Rplanet= "<< planet_radius << " Rearth." << endl;
      }
      if (in_c == 2) {
        cout << "Reading dust parameters... " << endl;
        A_dust = stod(line.substr(11,17));
        B_dust = stod(line.substr(22,26));
        alpha_dust = stod(line.substr(38,42));
        rho_dust = stod(line.substr(53,57));
        s_0 = stod(line.substr(65,70)) * 1.0e-4;
        cout << "rho= " << rho_dust << "g/cm3, initial size= " << s_0 << " micron." << endl;
        cout << "Sublimation coefficient= " << alpha_dust << endl;
        cout << "Clausius-Claperyon parameters, A= " << A_dust << ", B= " << B_dust << endl;
        
      }
      if (in_c == 3) {
        rmin = stod(line.substr(14,18));
        //cout << rmin << endl;
        rmax = stod(line.substr(27,31));
        //cout << rmax << endl;
        tmin = stod(line.substr(40,44));
        //cout << tmin << endl;
        tmax = stod(line.substr(52,56));
         //cout << tmax << endl;
        pmin = stod(line.substr(65,69)) * -1.0;
         //cout << pmin << endl;
        pmax = stod(line.substr(78,82));
         //cout << pmax << endl;

      }
      if (in_c == 4) {
        outflow = stoi(line.substr(14));
        //cout << outflow << endl;
        mdot_read = stod(line.substr(23,26));
        //cout << mdot_read << endl;
        tau_type = stoi(line.substr(36));
        //cout << tau_type << endl;
        if (outflow==1)  {
          cout << "Outflow is radially outwards from the whole planetary surface." << endl;
        } else {
          cout <<  "Outflow is from the planet's dayside." << endl;
        }
        cout << "The planetary mass loss rate is " << mdot_read << " Mearth/Gyr." << endl;
        if (tau_type ==1) {
          cout << "Optical depth of dust is traced." << endl;
          bool tau_constant = false;
        } else {
          cout << "Optical depth of dust is constant and thin (0.1)." << endl;
          bool tau_constant = true;
        }
      }
      if (in_c == 5) {
        major_timestep = stod(line.substr(10,14));
        cout << major_timestep << endl;
        no_orbits = stod(line.substr(25,28));
        cout << no_orbits << endl;
        nparticles = stoi(line.substr(38,41));
        cout << nparticles << endl;
      }
      
      in_c = in_c + 1;
      
    }
    input.close();
  }


  long int total_particles = nparticles; //initial number of particles to start simulation with
  double t_common = major_timestep;
  double big_step = major_timestep; //big time step (in terms of planetary orbits)
  double end_t = no_orbits; // end time of simulation
  double total_t = 0.0; // total time that has passed, so 0 in the beginning
  long int current_particles = 0; // number of current particles in simulation

  if (tau_constant == false) {
  build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min); //grid for optical depth calculations
  cout << "Built grid for optical depth tracing." << endl;
  }
  add_particles(particles, current_particles, total_particles, 0.0); // call function that adds particles to simulation (in add_rm.cpp file)
  
  rm_particles(particles);
  //Solve particles is the main routine of the program. It is in particles.cpp.
  solve_particles(0.00, end_t, particles, total_particles,current_particles);
  return 0;

}
