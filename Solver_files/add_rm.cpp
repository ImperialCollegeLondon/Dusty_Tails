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
ofstream ofile("kic_1255b_035_spherical.bin", ios::out | ios::binary);

ofstream ray_tracer("ray_tracer_test_300cells.bin", ios::out | ios::binary);


//define grid limits for ray tracing calculation
//atm grid is uniform
//anything with "a" in defines the A grid, anything with "b" defines B grid
double r_a [301];
double r_b [300];
double theta_a [301];
double theta_b [300];
double phi_a [301];
double phi_b [300];

double n_cells = 300.0; //number of grid cells in each direction
double r_min = 0.0;
double r_max = 300.0;
double theta_min = 0;
double theta_max =  300.;
double phi_min = 0.0;
double phi_max = 300.;

//spacing of grid
double dr = (r_max - r_min)/ n_cells;
double dtheta = (theta_max - theta_min ) / n_cells;
double dphi = ( phi_max - phi_min) / n_cells;

//define 3d arrays to store extinctions and optical depths at each grid cell
double extinction [300][300][300] = {};
double optical_depth [300][300][300] = {};

//limits of the particle distribution: this needs to be checked if using a
//different planet or different initial conditions for particles
double d_r_min = 0.92;
double d_r_max = 1.18;
double d_t_min = 1.4;
double d_t_max = 1.8;
double d_p_min = -0.6;
double d_p_max = 0.05;

//spacing of grid cells in scale of particle distribution
double d_dr = (d_r_max - d_r_min)/ n_cells;
double d_dtheta = (d_t_max - d_t_min ) / n_cells;
double d_dphi = ( d_p_max - d_p_min) / n_cells;


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



                 grain.p_size = 0.35e-4; //initial grain size
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


void solve_particles(double total_t, double end_t, vector <Particle>& particles, \
                     long int total_particles, double t_common, double big_step, \
                    long int current_particles){

  double plot_time = 0.01;
  vector <double> updated_vector(8);

  while (total_t < end_t) {

    for( Particle& p : particles) {

        double no_particles = particles.size();

        ofile.write((char*) &total_t, sizeof(double));

        //ofile.write((char*) &no_particles, sizeof(long int));
        //cout << "particles " << particles.size() << endl;
        ofile.write((char*) &p.id, sizeof(long int));

        ofile.write((char*) &p.position[0], sizeof(double));

        ofile.write((char*) &p.position[1], sizeof(double));

        ofile.write((char*) &p.position[2], sizeof(double));

        ofile.write((char*) &p.p_size, sizeof(double));

        ofile.write((char*) &p.p_mass, sizeof(double));

        updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
        p.velocity[0], p.velocity[1], p.velocity[2], p.p_size}, total_t, t_common, p.h_updated);
        p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
        p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
        p.p_size = updated_vector[6];
        p.p_mass = dust_mass(p.p_size);
        p.h_updated = updated_vector[7];

    }

    rm_particles(particles);

    cout << " at " << total_t << endl;

    if ( total_t >= 0.49 ) {

      cout << "now at grid builder " << endl;

      build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min);
      calculation_ext(particles, extinction);
      optical_depth_test(extinction, optical_depth);

      for (unsigned int j = 0; j <300; j++){
                for (unsigned int k = 0; k < 300; k++){
                    double theta_local;
                    double phi_local;

                    theta_local = theta_reverse(theta_a[j]);
                    phi_local = phi_reverse(phi_a[k]);
                    ray_tracer.write((char*) &theta_local, sizeof(double));
                    ray_tracer.write((char*) &phi_local, sizeof(double));
                    ray_tracer.write((char*) &extinction [299][j][k], sizeof(double));
                    ray_tracer.write((char*) &optical_depth [299][j][k], sizeof(double));


            }
        }

    }

    if (total_t > plot_time) {
      current_particles = total_particles;
      total_particles = total_particles + 300;
      add_particles(particles, current_particles, total_particles, total_t);
      plot_time = plot_time + 0.01;
     }


    total_t = total_t + big_step;
    t_common = t_common + big_step;


  }
}
