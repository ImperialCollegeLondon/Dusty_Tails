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
#include <cstring>
#include "opacities.h"



using namespace std;
unsigned seed = 123;
mt19937 generator (seed); //set the seed for the distrubutions of particle positions

//open files to write data for python plotting
ofstream ofile("./data/output.bin", ios::out | ios::binary);
//ofstream ray_tracer("./data/grid_test.bin", ios::out | ios::binary);


//spacing of grid cells in scale of particle distribution
double d_dr = (d_r_max - d_r_min)/ r_cells_d;
double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
double d_dphi = ( d_p_max - d_p_min) / p_cells_d;

int no_particles[r_cells][t_cells][p_cells] = {};

vector < vector < tk:: spline > > s_phi;

//Use something of the sort below to have particles starting on the day side of the planet
uniform_real_distribution<double> uniform_phi_short(0.375, 0.625);
uniform_real_distribution<double> uniform_theta_short(0.2, 0.8);
//Use distrubutions below to have particles coming out of the whole planetary surface
uniform_real_distribution<double> uniform_phi(0.0, 1.0);
uniform_real_distribution<double> uniform_theta(0.0, 1.0);


//function add particles adds more particles to the simulation at a given time
//as arguments it takes the vector of particles, the current number of particles,
//the total of particles we want to get and the current time in the simulation
void add_particles(vector <Particle> &particles, long int current,
                   long int total, double time){

    Particle grain;
    double kappa_init =  opac.stellar_abs(s_0) 
                     + opac.stellar_scat(s_0);
    for ( unsigned long int i = current; i < total; i++){
        double phi,theta;
        grain.id = i+1; //number ID of particle
        
        //generate particles initial position
        if (outflow==1) {
        //spherical
        phi = 2.0*PI * uniform_phi(generator); //phi coordinate
        theta = acos(1- 2.0* uniform_theta(generator)); //theta coordinate
        } else {
        //dayside
        phi = 2.0*PI * uniform_phi_short(generator); //phi coordinate
        theta = acos(1- 2.0* uniform_theta_short(generator)); //theta coordinate
        }

        //initial position of particle in cartesian
        grain.position = {r_start*sin(theta)*cos(phi) + planet_x, \
                        r_start*sin(theta)*sin(phi), \
                        r_start*cos(theta)};

        //velocity of particle in cartesian
        grain.velocity = {v_esc*sin(theta)*cos(phi), \
                        v_esc*sin(theta)*sin(phi), \
                        v_esc*cos(theta)};

        grain.p_size = s_0; //initial grain size

        //Initial particle optical depth
        if (tau_constant == true) {
            grain.p_tau = 0.1;
        } else  {
            grain.p_tau = 0.001;
        }
        
        grain.h_updated = 1.0e-6; //initial time step for numerical integrator
        grain.p_mass = dust_mass(grain.p_size); //initial particle mass
        grain.p_opacity = kappa_init; //initial rad pressure efficiency
        //initial particle temperature
        grain.p_temp = brent( grain.p_size, grain.position[0], 
                                grain.position[1], grain.position[2], grain.p_tau); 
        //cout << grain.p_temp << endl;
        particles.push_back(grain); //add particle to the vector of particles
                 

    }

}

//The function below removes particles that have become too small 
//for it to be worth tracking
//Argument is the vector of particles
void rm_particles(vector <Particle>& particles){
    for (unsigned long int i = 0; i < particles.size(); i++){
      if (particles[i].id == 14) {
          cout << particles[i].p_size << endl;
      }
      if (particles[i].p_size < 0.01e-4) {
        
        particles.erase(particles.begin() + i);
        if (i!=0) {
        i--; }
      }
      if (isnan(particles[i].p_size)  ) {
        cout << "Size of particle is NaN, something has gone wrong " << endl;
        particles.erase(particles.begin() + i);
        i--;
      }
    }
}

// solve_particles takes as arguments the total time of the simulation (total_t)
// the end time of the current iteration (end_t), the vector of no_particles
// the total number of particles to be achieved (total_particles)
//the current number of particles in the iteration (current_particles)

void solve_particles(double total_t, double end_t, vector <Particle>& particles, \
                     long int total_particles, long int current_particles){
  double previous_t = 0.00;
  double plot_time = 0.01;
  double t_next = 0.00; 
  //vector which will take updated values of positons, velocitites, size and optimal time step for particle
  vector <double> updated_vector(8); 

  //variables related to the optical depth interpolation
  vector < vector <double> > s_phi_values;
  vector < tk:: spline > s_theta;
  vector <double> s_theta_values;
  vector < vector < vector <double> > >  tau;
  vector <double> radii_v, thetas_v, phis_v;

  //Intepolator needs grid in vector format.
  if (tau_constant == false ) {
    radii_v  = r_grid_to_vector(r_a);
    thetas_v = t_grid_to_vector(theta_a);
    phis_v   = p_grid_to_vector(phi_a);
   }


  while (t_next < end_t) {
    double t_global_min = 0.01;
    cout << "At orbit " << t_next << endl;
    if ((tau_constant == false) && (t_next > 1.0)) {
        memset(extinction, 0.0, sizeof(extinction));
        memset(optical_depth, 0.0, sizeof(optical_depth));
        calculation_ext(particles, extinction, t_global_min);
        optical_depth_calc(extinction, optical_depth);
        tau.clear();
        tau = tau_to_vector(optical_depth);
        s_phi.clear();
        s_phi = splines_phi(tau, radii_v, thetas_v, phis_v);

        t_global_min = 1.0e-3;
    } 
  //do bigger time steps in first orbit and smaller ones afterwards
  //to establish a primary population of dust particles
  t_next = t_global_min + previous_t;
  if (t_next >= plot_time) {
           t_next = plot_time;
  }
  int tau_thick = 0; //counter for particles with a thick optical depth
  
  for( Particle& p : particles) {
    double no_particles = particles.size();
    vector <double> pos_spherical;
    vector <double> pos_scaled;
    if ((tau_constant == true) || (previous_t <= 1.0) ) {
            p.p_tau = 0.001;
    } else {
            pos_spherical = pos_to_spherical(p.position[0], p.position[1], p.position[2]);
            pos_scaled = grid_scaling(pos_spherical); //position of particle in optical depth grid
            if (pos_scaled[0] == -1.0) {
                //if particle sits outside the grid the optical depth is set to 1e-3
                p.p_tau = 0.001;
            } else {
                //if particle sits inside the grid the optical depth of it is 
                //calculated via a trilinear spline interpolation
                s_phi_values = phi_spline_result(s_phi , radii_v, thetas_v, phis_v, pos_scaled[2]);
                s_theta =  splines_theta ( s_phi_values, radii_v, thetas_v);
                s_theta_values = theta_spline_result( s_theta, radii_v, thetas_v, phis_v , pos_scaled[1]);
                p.p_tau =  tau_p (s_theta_values, radii_v, thetas_v, phis_v, pos_scaled[0]);
                
                if (p.p_tau <=0.001) {
                    p.p_tau = 0.001;}
    
                if (p.p_tau > 1.0){
                    tau_thick = tau_thick + 1;}
            }
            s_phi_values.clear();
            s_theta.clear();
            s_theta_values.clear();
        }
        //update all relevant particle components
        //function RK solver calls the numerical solver to obtain the new particle
        //positions, velocities, size (and mass - consequence of the size change)
        updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
        p.velocity[0], p.velocity[1], p.velocity[2], p.p_size, p.p_tau, p.p_temp}, previous_t, t_next, p.h_updated);
        p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
        p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
        p.p_size = updated_vector[6];
        p.p_opacity = opac.stellar_abs(p.p_size) + opac.stellar_scat(p.p_size);
        p.p_mass = dust_mass(p.p_size); 
        p.h_updated = updated_vector[8];    
        p.p_temp = brent( p.p_size, p.position[0], 
                                p.position[1], p.position[2], p.p_tau); 
        //cout << p.id << endl;
    }
    cout << "Number of optically thick particles " << tau_thick << endl;
    rm_particles(particles); //removes particles that are too small
    
    //to obtain optical depth plot when optical depth is constant
    /*
    if (t_next > (end_t-t_global_min)){
            calculation_ext(particles, extinction, t_global_min);
            optical_depth_calc(extinction, optical_depth);

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
        */

    if (abs(t_next-plot_time) < 1.0e-8) {  
      for( Particle& p : particles) {
        //write file with particle data
        ofile.write((char*) &t_next, sizeof(double));
        ofile.write((char*) &p.id, sizeof(long int));
        ofile.write((char*) &p.position[0], sizeof(double));
        ofile.write((char*) &p.position[1], sizeof(double));
        ofile.write((char*) &p.position[2], sizeof(double));
        ofile.write((char*) &p.p_size, sizeof(double));
        ofile.write((char*) &p.p_mass, sizeof(double));
        ofile.write((char*) &p.p_tau, sizeof(double));
        ofile.write((char*) &p.p_temp, sizeof(double));
        ofile.write((char*) &p.p_opacity, sizeof(double));
      }
      current_particles = total_particles;
      total_particles = total_particles + 1000;

      //add particles every 100th of an orbit
      add_particles(particles, current_particles, total_particles, t_next);
      plot_time = plot_time + 0.01;
     }
    previous_t = t_next;
    
  }
}

//*******************************************************************************
//*******************************************************************************
//interpolation routine for the optical depth - functions are written in spline.h

double cubicInterpolate ( vector <double> p, vector <double> s1, double x) {
    tk::spline s(s1,p, tk::spline::cspline);
    return s(x);
    
}

double bicubicInterpolate (vector< vector <double> > p, vector <double> s1, 
                                    vector <double> s2, double x, double y) {
    vector <double> arr(s1.size());
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = cubicInterpolate(p[i], s2, y);
     }
    return cubicInterpolate(arr, s1, x);
}

double tricubicInterpolate (vector <vector< vector <double> > > p,
                     vector <double> s1, vector <double> s2, vector <double> s3, 
                     double x, double y, double z) {

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



vector < vector < tk:: spline > > splines_phi (vector <vector< vector <double> > > taus, 
                    vector <double> radii, vector <double> thetas, vector <double> phis) {

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
vector <vector <double> > phi_spline_result(vector < vector < tk:: spline >> splines , 
                            vector<double> radii, \
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

vector < tk:: spline > splines_theta( vector <vector <double>> phi_splines_p, 
                                      vector<double> radii, vector<double> thetas) {
     
    vector < tk::spline > r_splines;
    for (int i =0; i < radii.size(); i++){
        tk::spline s(thetas, phi_splines_p[i], tk::spline::cspline);
        r_splines.push_back(s);
    }

    return r_splines;
}
vector <double> theta_spline_result( vector <tk:: spline> splines, 
                vector<double> radii, vector<double> thetas, 
                vector<double> phis , double theta) {

    vector <double>  theta_spline_values;
    for (int i = 0; i<radii.size(); i++){
            theta_spline_values.push_back(splines[i](theta));}

    return theta_spline_values;
    } 

 double tau_p (vector<double> theta_splines_p, vector<double> radii, 
               vector<double> thetas, vector<double> phis, double radius){

    tk:: spline s(radii, theta_splines_p);
    return s(radius);
}
//*******************************************************************************
//*******************************************************************************
