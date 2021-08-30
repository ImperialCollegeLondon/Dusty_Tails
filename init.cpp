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
double r_a [r_cells + 1];
double r_b [r_cells];
double theta_a [t_cells + 1];
double theta_b [t_cells];
double phi_a [p_cells + 1];
double phi_b [p_cells];


double r_min = 0.0;
double r_max = r_cells_d;
double theta_min = 0.0;
double theta_max =  t_cells_d;
double phi_min = 0.0;
double phi_max = p_cells_d;

double dr = (r_max - r_min)/ r_cells_d;
double dtheta = (theta_max - theta_min ) / t_cells_d;
double dphi = ( phi_max - phi_min) / p_cells_d;

  double extinction [r_cells][t_cells][p_cells];
  double optical_depth [r_cells][t_cells][p_cells];


int main() {

  long int total_particles = 1000; //initial number of particles to start simulation with
  double t_common = 0.01;
  double big_step = 0.01; //big time step (in terms of planetary orbits)
  double end_t = 5.0; // end time of simulation
  double total_t = 0.0; // total time that has passed, so 0 in the beginning


  long int current_particles = 0; // number of current particles in simulation

//define grid limits for ray tracing calculation
//atm grid is uniform
//anything with "a" in defines the A grid, anything with "b" defines B grid
 
  build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min); //grid for optical depth calculations
  add_particles(particles, current_particles, total_particles, 0.0); // call function that adds particles to simulation (in add_rm.cpp file)
  cout << "built grids and added particles" << endl;
  rm_particles(particles);
  cout << "rm particles fine" << endl;
  calculation_ext(particles, extinction);
  optical_depth_calc(extinction, optical_depth);
  solve_particles(0.00, end_t, particles, total_particles,current_particles);
  return 0;

}
