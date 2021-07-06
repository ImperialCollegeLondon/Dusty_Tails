#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include "spline.h"
#include <random>
#include <chrono>
#include <stdlib.h>

using namespace std;


unsigned seed = 123;
mt19937 generator (seed); //set the seed for the distrubutions of particle positions

//Use something of the sort below to have particles starting on the day side of the planet
uniform_real_distribution<double> uniform_phi(0.375, 0.625);
uniform_real_distribution<double> uniform_theta(0.2, 0.8);

//Use distrubutions below to have particles coming out of the whole planetary surface
//uniform_real_distribution<double> uniform_phi(0.0, 1.0);
//uniform_real_distribution<double> uniform_theta(0.0, 1.0);

//open files to write data for python plotting
ofstream ofile("./data/KIC1255b_5orb_035_tau_01micro_1node.bin", ios::out | ios::binary);

ofstream ray_tracer("./data/ray_KIC1255b_5orb_035_tau_01micro_1node.bin", ios::out | ios::binary);


//define 3d arrays to store extinctions and optical depths at each grid cell
//double extinction [r_cells][t_cells][p_cells] = {};
//double optical_depth [r_cells][t_cells][p_cells] = {};
int no_particles[r_cells][t_cells][p_cells] = {};

//limits of the particle distribution: this needs to be checked if using a
//different planet or different initial conditions for particles
//this is obviously not very sustainable
double d_r_min = 0.9;
double d_r_max = 1.1;
double d_t_min = 1.55;
double d_t_max = 1.60;
double d_p_min = -0.25;
double d_p_max = 0.0;

//spacing of grid cells in scale of particle distribution
double d_dr = (d_r_max - d_r_min)/ r_cells_d;
double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
double d_dphi = ( d_p_max - d_p_min) / p_cells_d;

vector < vector < tk:: spline > > s_phi;


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

                 grain.pos_spherical = pos_to_spherical(grain.position[0], grain.position[1], grain.position[2]);

                 grain.v_spherical = vel_to_spherical(grain.velocity[0], grain.velocity[1], grain.velocity[2]);

                 grain.p_size = 0.35e-4; //initial grain size
                 grain.p_tau = 0.0;
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
  vector < vector <double> > s_phi_values;
  vector < tk:: spline > s_theta;
  vector <double> s_theta_values;
  vector < vector < vector <double> > >  tau;
  vector <double> radii_v, thetas_v, phis_v;

  radii_v = r_grid_to_vector(r_a);
  thetas_v = t_grid_to_vector(theta_a);
  phis_v = p_grid_to_vector(phi_a);

  while (total_t < end_t) {
    double extinction [r_cells][t_cells][p_cells] = {};
    double optical_depth [r_cells][t_cells][p_cells] = {};
      //calculation_ext is in ray_tracer.cpp - calculates the extinction at each grid cell
  calculation_ext(particles, extinction);
  //optical_depth_test is in ray_tracer.cpp - calculates the optical depth in each grid cell, dependent on the extinction distribution
  optical_depth_calc(extinction, optical_depth);
  
  tau.clear();
  tau = tau_to_vector(optical_depth);
  s_phi.clear();
  s_phi = splines_phi( tau, radii_v, thetas_v, phis_v);

  cout << "orbit: " << total_t << endl;

  for( Particle& p : particles) {

        double no_particles = particles.size();
        vector <double> pos_scaled;
        pos_scaled = grid_scaling(p.pos_spherical);
        
        if (pos_scaled[0] == -1.0) {
            p.p_tau = 0.0;
        } else {
        s_phi_values = phi_spline_result(s_phi , radii_v, thetas_v, phis_v, pos_scaled[2]);
        s_theta =  splines_theta ( s_phi_values, radii_v, thetas_v);
        s_theta_values = theta_spline_result( s_theta, radii_v, thetas_v, phis_v , pos_scaled[1]);
        p.p_tau =  tau_p (s_theta_values, radii_v, thetas_v, phis_v, pos_scaled[0]);
        if (p.p_tau <=1.e-33) {
          p.p_tau = 0.0;
        }
        if (isnan(p.p_tau) ) {
            p.p_tau = 0.0;
        }
        }
        
        //Lines to write binary file
        ofile.write((char*) &total_t, sizeof(double));
        ofile.write((char*) &p.id, sizeof(long int));
        ofile.write((char*) &p.position[0], sizeof(double));
        ofile.write((char*) &p.position[1], sizeof(double));
        ofile.write((char*) &p.position[2], sizeof(double));
        ofile.write((char*) &p.p_size, sizeof(double));
        ofile.write((char*) &p.p_mass, sizeof(double));
        ofile.write((char*) &p.p_tau, sizeof(double));

        s_phi_values.clear();
        s_theta.clear();
        s_theta_values.clear();

        //updated vector is new positions, velocities, size and optimal small step size for particle
        // RK_solver function is in "solver_new_err.cpp"
        updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
        p.velocity[0], p.velocity[1], p.velocity[2], p.p_size, p.p_tau}, total_t, t_common, p.h_updated);
        p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
        p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
        p.p_size = updated_vector[6];
        p.p_mass = dust_mass(p.p_size); //dust_mass is in "microphysics.cpp"
        p.h_updated = updated_vector[7];
        p.pos_spherical = pos_to_spherical(p.position[0], p.position[1], p.position[2]);
        p.v_spherical = vel_to_spherical(p.velocity[0], p.velocity[1], p.velocity[2]);

    }

    rm_particles(particles); //removes particles that are too small

    //if condition below is just for a test of the ray tracer at a given time
    if ( total_t >= 4.99 ) {
      //calculation_ext is in ray_tracer.cpp - calculates the extinction at each grid cell
      //calculation_ext(particles, extinction, no_particles);
      //optical_depth_test is in ray_tracer.cpp - calculates the optical depth in each grid cell, dependent on the extinction distribution
      //optical_depth_test(extinction, optical_depth);
      //tau = tau_to_vector(optical_depth);


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


    //when plot time is reached the following condition adds new no_particles
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

double cubicInterpolate ( vector <double> p, vector <double> s1, double x) {
    tk::spline s(s1,p, tk::spline::cspline);
    return s(x);
    
}

double bicubicInterpolate (vector< vector <double> > p, vector <double> s1, vector <double> s2, double x, double y) {
    vector <double> arr(s1.size());
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = cubicInterpolate(p[i], s2, y);
     }
    return cubicInterpolate(arr, s1, x);
}

double tricubicInterpolate (vector <vector< vector <double> > > p, vector <double> s1, vector <double> s2, vector <double> s3, double x, double y, double z) {
    vector <double> arr(s1.size());
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = bicubicInterpolate(p[i], s2, s3, y, z);
     }
    double tri = cubicInterpolate(arr, s1, x);
    if (tri < 1.0e-20) {
        return 0.0;
    
    } else {
        return tri;
    }
}

vector < vector < tk:: spline > > splines_phi (vector <vector< vector <double> > > taus, vector <double> radii, vector <double> thetas, vector <double> phis) {
    vector < vector < tk::spline >> phi_splines;
    tk::spline s; 
    for (int i = 0; i<radii.size(); i++){
        phi_splines.push_back({});
        for ( int j=0; j<thetas.size(); j++) {
            tk::spline s(phis, taus[i][j], tk::spline::cspline);
            phi_splines[i].push_back(s);
        }
    }

    return phi_splines;
}
vector <vector <double> > phi_spline_result(vector < vector < tk:: spline >> splines , vector<double> radii, \
                            vector <double> thetas, vector<double> phis, double phi) {
    vector <vector <double> > phi_spline_values;
    for (int i = 0; i<radii.size(); i++){
        phi_spline_values.push_back({});
        for ( int j=0; j<thetas.size(); j++) {
            phi_spline_values[i].push_back(splines[i][j](phi));
        }
    } 
    return phi_spline_values;
}

vector < tk:: spline > splines_theta( vector <vector <double>> phi_splines_p, vector<double> radii, vector<double> thetas) {
     
    vector < tk::spline > r_splines;
    for (int i =0; i < radii.size(); i++){
        tk::spline s(thetas, phi_splines_p[i], tk::spline::cspline);
        r_splines.push_back(s);
    }

    return r_splines;
}
vector <double> theta_spline_result( vector <tk:: spline> splines, vector<double> radii, vector<double> thetas, vector<double> phis , double theta) {
    vector <double>  theta_spline_values;
    for (int i = 0; i<radii.size(); i++){
            theta_spline_values.push_back(splines[i](theta));
       
        }
    return theta_spline_values;
    } 

 double tau_p (vector<double> theta_splines_p, vector<double> radii, vector<double> thetas, vector<double> phis, double radius){
    tk:: spline s(radii, theta_splines_p);
    
    return s(radius);
}

