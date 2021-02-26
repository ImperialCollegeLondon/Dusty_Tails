#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include <random>
#include <chrono>

using namespace std;

unsigned seed = 123;
mt19937 generator (seed); //set the seed for the distrubutions of particle positions

//Use something of the sort below to have particles starting on the day side of the planet
//uniform_real_distribution<double> uniform_phi(0.625, 0.875);
//uniform_real_distribution<double> uniform_theta(0.2, 0.8);

//Use distrubutions below to have particles coming out of the whole planetary surface
uniform_real_distribution<double> uniform_phi(0.0, 1.0);
uniform_real_distribution<double> uniform_theta(0.0, 1.0);

//open files to write data for python plotting
ofstream ofile("KIC1255b_test.bin", ios::out | ios::binary);

ofstream ray_tracer("ray_tracer_test.bin", ios::out | ios::binary);


//define grid limits for ray tracing calculation
//atm grid is uniform
//anything with "a" in defines the A grid, anything with "b" defines B grid
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

//spacing of grid
double dr = (r_max - r_min)/ r_cells_d;
double dtheta = (theta_max - theta_min ) / t_cells_d;
double dphi = ( phi_max - phi_min) / p_cells_d;

//define 3d arrays to store extinctions and optical depths at each grid cell
double extinction [r_cells][t_cells][p_cells] = {};
double optical_depth [r_cells][t_cells][p_cells] = {};
int no_particles[r_cells][t_cells][p_cells] = {};

//limits of the particle distribution: this needs to be checked if using a
//different planet or different initial conditions for particles
//this is obviously not very sustainable
double d_r_min = 0.9;
double d_r_max = 1.1;
double d_t_min = 1.55;
double d_t_max = 1.60;
double d_p_min = -0.3;
double d_p_max = 0.0;

//spacing of grid cells in scale of particle distribution
double d_dr = (d_r_max - d_r_min)/ r_cells_d;
double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
double d_dphi = ( d_p_max - d_p_min) / p_cells_d;


//function add particles adds more particles to the simulation at a given time
//as arguments it takes the vector of particles, the current number of particles,
//the total of particles we want to get and the current time in the simulation
void add_particles(vector <Particle> &particles, long int current,
                   long int total, double time){

        double v_esc; //escape velocity
        Particle grain;

        v_esc = (pow((2.0*G *0.05*Mearth)/(0.38*Rearth), 0.5)) * (T/a);
        //escape velocity m/s

        for ( unsigned long int i = current; i < total; i++){

                 //cout << "adding particle " << i << endl;

                 grain.id = i+1; //number ID of particle
                 double phi = 2.0*PI * uniform_phi(generator); //phi coordinate

                 double theta = acos(1- 2.0* uniform_theta(generator)); //theta coordinate

                 //position of particle in cartesian
                 grain.position = {r_start*sin(theta)*cos(phi) + planet_x, \
                                   r_start*sin(theta)*sin(phi), \
                                   r_start*cos(theta)};

                //velocity of particle in cartesian
                 grain.velocity = {v_esc*sin(theta)*cos(phi), \
                                   v_esc*sin(theta)*sin(phi), \
                                   v_esc*cos(theta)};


                 grain.p_size = 0.40e-4; //initial grain size
                 grain.p_tau = tau; //using a constant optical depth for now defined in constants.h
                 grain.p_density = rho_d; //bulk density
                 grain.h_updated = 0.001; //initial time step for numerical integrator
                 grain.p_mass = dust_mass(grain.p_size); //initial particle mass

                 particles.push_back(grain); //add particle to the vector of particles

               }

}

//The function below removes particles that have become too small for it to be worth of being considered
//Argument is the vector of particles
void rm_particles(vector <Particle>& particles){
    for (unsigned long int i = 0; i < particles.size(); i++){
      if (particles[i].p_size < 0.1e-4) {
        particles.erase(particles.begin() + i);
        i--;
      }
    }
}

// solve_particles takes as arguments the total time of the simulation (total_t)
// the end time of the current iteration (end_t), the vector of no_particles
// the total number of particles to be achieved (total_particles)
//the current number of particles in the iteration (current_particles)
//t_common and big_step are the big common time step
void solve_particles(double total_t, double end_t, vector <Particle>& particles, \
                     long int total_particles, double t_common, double big_step, \
                    long int current_particles){

  double plot_time = 0.01; //time when to output values for plotting
  vector <double> updated_vector(8); //vector which will take updated values of positons, velocitites, size and optimal time step for particle
  int counter_ps = 0;
  double optical_depth_avg = 0.0;
  double optical_depth_sum = 0.0;
  double particles_in_box = 0.;
  double no_particles_avg = 0.0;

  cout << "dr " << d_dr << endl;
  cout << "dtheta " << d_dtheta << endl;
  cout << "dphi " << d_dphi << endl;

  while (total_t < end_t) {

    cout << "orbit: " << total_t << endl;

    for( Particle& p : particles) {

        double no_particles = particles.size();

        //Lines to write binary file
        ofile.write((char*) &total_t, sizeof(double));
        ofile.write((char*) &p.id, sizeof(long int));
        ofile.write((char*) &p.position[0], sizeof(double));
        ofile.write((char*) &p.position[1], sizeof(double));
        ofile.write((char*) &p.position[2], sizeof(double));
        ofile.write((char*) &p.p_size, sizeof(double));
        ofile.write((char*) &p.p_mass, sizeof(double));

        //updated vector is new positions, velocities, size and optimal small step size for particle
        // RK_solver function is in "solver_new_err.cpp"
        updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
        p.velocity[0], p.velocity[1], p.velocity[2], p.p_size}, total_t, t_common, p.h_updated);
        p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
        p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
        p.p_size = updated_vector[6];
        p.p_mass = dust_mass(p.p_size); //dust_mass is in "microphysics.cpp"
        p.h_updated = updated_vector[7];

    }

    rm_particles(particles); //removes particles that are too small

    //if condition below is just for a test of the ray tracer at a given time

    if ( total_t >= 0.49 ) {

      cout << "now at grid builder " << endl;
      //build_grids is in ray_tracer.cpp - as the name says it builds the grid over the star for the ray tracing calculations
      build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min);
      //calculation_ext is in ray_tracer.cpp - calculates the extinction at each grid cell
      calculation_ext(particles, extinction, no_particles);
      //optical_depth_test is in ray_tracer.cpp - calculates the optical depth in each grid cell, dependent on the extinction distribution
      optical_depth_test(extinction, optical_depth);

      for (unsigned int l = 0; l <r_cells; l++){
        for (unsigned int m = 0; m <t_cells; m++){
          for (unsigned int n = 0; n <p_cells; n++){
            optical_depth_sum = optical_depth_sum + optical_depth[l][m][n];
            //cout << "at no particles writer" << endl;
            if (no_particles[l][m][n] != 0) {
              if (no_particles[l][m][n] == 1 ) {
                cout << "i " << l << endl;
                cout << "j " << m << endl;
                cout << "k " << n << endl;
                cout << "no particles in cell " << no_particles[l][m][n] << endl;
                cout << "optical depth " << optical_depth[l][m][n] << endl; }

              particles_in_box = particles_in_box + no_particles[l][m][n]; }
              else {
                counter_ps = counter_ps + 1; }
              }
              }
        }

      //loop below write file for plotting
      for (unsigned int j = 0; j <t_cells; j++){
                for (unsigned int k = 0; k < p_cells; k++){
                    double theta_local;
                    double phi_local;

                    theta_local = theta_reverse(theta_a[j]);
                    phi_local = phi_reverse(phi_a[k]);
                    ray_tracer.write((char*) &theta_local, sizeof(double));
                    ray_tracer.write((char*) &phi_local, sizeof(double));
                    ray_tracer.write((char*) &extinction [r_cells-1][j][k], sizeof(double));
                    ray_tracer.write((char*) &optical_depth [r_cells-1][j][k], sizeof(double));


            }
        }
      }


    //when plot time is reached the following condition adds 100 new no_particles
    //should change this 100 to a variable
    if (total_t > plot_time) {
      current_particles = total_particles;
      total_particles = total_particles + 1000;
      //add particles explained above in this file
      add_particles(particles, current_particles, total_particles, total_t);
      plot_time = plot_time + 0.01;
     }


    total_t = total_t + big_step;
    t_common = t_common + big_step;


  }
}
