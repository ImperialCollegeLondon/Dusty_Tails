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
#include <omp.h>


#include <chrono>
using namespace std::chrono;


using namespace std;
unsigned seed = 123;
mt19937 generator (seed); //set the seed for the distrubutions of particle positions
long int counter = 0;

//open files to write data for python plotting
ofstream output("./data/output.bin", ios::out | ios::binary);
ofstream output_tau("./data/tau_data.bin", ios::out | ios::binary);
//ofstream output_tau_test("./data/tau_test.bin", ios::out | ios::binary);
//ofstream ray_tracer("./data/grid_test.bin", ios::out | ios::binary);
ofstream output_final("./data/output_final.bin", ios::out | ios::binary);

//spacing of grid cells in scale of particle distribution
double d_dr = (d_r_max - d_r_min)/ r_cells_d;
double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
double d_dphi = ( d_p_max - d_p_min) / p_cells_d;

vector < vector < tk:: spline > > s_phi;

//Use something of the sort below to have particles starting on the day side of the planet
uniform_real_distribution<double> uniform_phi_short(0.375, 0.625);
uniform_real_distribution<double> uniform_theta_short(0.2, 0.8);
//Use distrubutions below to have particles coming out of the whole planetary surface
uniform_real_distribution<double> uniform_phi(0.0, 1.0);
uniform_real_distribution<double> uniform_theta(0.0, 1.0);

std::normal_distribution<double> ndist{2.00,1.00};

vector <dust_read> read_data(){
  std::fstream output;
  output.open("input.bin", std::fstream::in | std::fstream::binary);
  output.seekg(0, ios::end);
  int size=output.tellg();
  output.seekg(0, ios::beg);
  long int total = size/sizeof(dust);
  vector <dust_read> dust_grains_in;
  for(int i = 0; i < total; i++){
       dust_grains_in.push_back(dust_read());
       output.read((char *) &dust_grains_in[i], sizeof(dust_read));
       

       //cout << dust_grains_out[i].id << endl;
    }
   output.close();
   return dust_grains_in;
}

//function add particles adds more particles to the simulation at a given time
//as arguments it takes the vector of particles, the current number of particles,
//the total of particles we want to get and the current time in the simulation
void add_particles(vector <Particle> &particles, long int &current_particles, long int &total_particles, double time){

    Particle grain;
    double opac_abs_init, opac_scat_init, kappa_planck_init, gsca_init;
    if (s_dist == 0) {
    double opac_abs_init = opac.stellar_abs(s_0);
    double opac_scat_init = opac.stellar_scat(s_0);
    double kappa_planck_init =  opac_abs_init + opac_scat_init;
    double gsca_init = opac.stellar_gsc(s_0);
    } 

    long int current, total;
     //flag for continuing dust
    vector <dust_read> particles_read;
    if (cont ==1) {
        //total number of particles is the size of the dust_grains_in vector
        cout << "Data from previous run is being read... " << endl;
        particles_read = read_data();
        
        //double previous_tstamp = 0.0;
        for (dust_read& p : particles_read) {
            grain.id = p.id;
            grain.position = {p.x_dust, p.y_dust, p.z_dust};
            grain.velocity = {p.vx_dust, p.vy_dust, p.vz_dust};
            grain.size = p.s_dust;
            grain.h_updated = p.h_dust;
            grain.mass = p.m_dust;
            grain.opac_planck = p.kappa_planck;
            grain.opac_abs = p.kappa_dust_abs;
            grain.opac_scat = p.kappa_dust_scat;
            grain.n_mini = p.nmini;
            grain.gsca = opac.stellar_gsc(grain.size);
            grain.temp_d = p.temp_dust;
            grain.pos_spherical = pos_to_spherical(grain.position[0], grain.position[1], grain.position[2]);
            grain.tau_d = p.tau_dust;
            particles.push_back(grain);
            
        }
    current = particles.size();
    cout << "Tail has " << current << " super-particles from previous run. " << endl;    
    }
    
    current = particles.size();
    int current_id = (time/major_timestep) *nparticles;
    total = current + nparticles;
    for ( unsigned long int i = current; i < total; i++){
        double phi,theta;
        //grain.id = i+1; //number ID of particle
        grain.id = current_id;
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
        grain.position = {(1.1*r_h)*sin(theta)*cos(phi) + planet_x, \
                        (1.1*r_h)*sin(theta)*sin(phi), \
                        (1.1*r_h)*cos(theta)};

        //velocity of particle in cartesian
        grain.velocity = {v_th*sin(theta)*cos(phi), \
                        v_th*sin(theta)*sin(phi), \
                        v_th*cos(theta)};

        if (s_dist == 0) {
        grain.size = s_0; //initial grain size
        grain.n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
        } else if (s_dist ==1) {
            double size_temp;
            size_temp = ndist(generator)* 1.0e-4;
            while (size_temp <= (0.01* 1.0e-4)) {
                size_temp = ndist(generator)* 1.0e-4;
            }
            grain.size = size_temp;
            grain.n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(grain.size, 3));
        }

        //Initial particle optical depth
        if (tau_constant == true) {
            grain.tau_d = 0.1;
        } else  {
            grain.tau_d = 0.001;
        }
        
        grain.h_updated = 1.0e-4; //initial time step for numerical integrator
        grain.mass = dust_mass(grain.size); //initial particle mass
        if (s_dist == 1) {
        double opac_abs_init = opac.stellar_abs(grain.size);
        double opac_scat_init = opac.stellar_scat(grain.size);
        double kappa_planck_init =  opac_abs_init + opac_scat_init;
        double gsca_init = opac.stellar_gsc(grain.size);
        } 
        grain.opac_planck = kappa_planck_init; //initial rad pressure efficiency
        //initial particle temperature
        grain.temp_d = brent( grain.size, grain.position[0], 
                                grain.position[1], grain.position[2], grain.tau_d); 
        grain.pos_spherical = pos_to_spherical(grain.position[0], grain.position[1], grain.position[2]);
        grain.opac_abs = opac_abs_init;
        grain.opac_scat = opac_scat_init;
        grain.gsca = gsca_init;
        
        particles.push_back(grain); //add particle to the vector of particles

        current_id = current_id +1;        
    }
    cont = 0;
    
    current_particles = current;
    total_particles = total;

}

//The function below removes particles that have become too small 
//for it to be worth tracking
//Argument is the vector of particles
bool _predicate(Particle& element) {
    
    return (element.size <= 1.0e-6); 
    }
void rm_particles(vector <Particle>& particles){

    particles.erase(std::remove_if(particles.begin(), 
    particles.end(), _predicate), particles.end());


    particles.shrink_to_fit();
    cout << "There are " << particles.size() << " super-particles in the tail." << endl;
}

// solve_particles takes as arguments the total time of the simulation (total_t)
// the end time of the current iteration (end_t), the vector of no_particles
// the total number of particles to be achieved (total_particles)
//the current number of particles in the iteration (current_particles)

void solve_particles(double total_t, double end_t, vector <Particle>& particles, \
                     long int total_particles, long int current_particles){
  double current_t = total_t;
  double plot_time = total_t;
  double t_next = total_t; 
  //vector which will take updated values of positons, velocitites, size and optimal time step for particle
  vector <double> updated_vector(8); 
  double d_dr = (d_r_max - d_r_min)/ r_cells_d;
  double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
  double d_dphi = ( d_p_max - d_p_min) / p_cells_d;


  //variables related to the optical depth interpolation
  vector < vector <double> > s_phi_values;
  vector < tk:: spline > s_theta;
  vector <double> s_theta_values;
  vector < vector < vector <double> > >  tau;
  vector <double> radii_v, thetas_v, phis_v;

  //Intepolator needs grid in vector format.
  
    radii_v  = r_grid_to_vector(r_a);
    thetas_v = t_grid_to_vector(theta_b);
    phis_v   = p_grid_to_vector(phi_b);


  while (current_t < end_t) {
    double t_global_min = major_timestep;
    auto start = high_resolution_clock::now();
    cout << "At time: " << current_t << endl;

    if (tau_constant == false) {
        memset(extinction, 0.0, sizeof(extinction));
        memset(optical_depth, 0.0, sizeof(optical_depth));
        calculation_ext(particles, extinction, t_global_min);
        optical_depth_calc(extinction, optical_depth);
         
        cout << "Obtained optical depth grid successfully." << endl;
        tau.clear();
        tau.shrink_to_fit();
        tau = tau_to_vector(optical_depth);
        s_phi.clear();
        s_phi.shrink_to_fit();
        s_phi = splines_phi(tau, radii_v, thetas_v, phis_v);

        t_global_min = 5.0e-3;
       
    } 


  t_next = t_global_min + current_t;
  std::vector<Particle>::iterator it;


  #pragma omp parallel for private( s_phi_values, s_theta, s_theta_values)
  for( Particle& p : particles) {

    vector <double> pos_spherical;
    vector <double> pos_scaled;
    if (tau_constant == true) {
            p.tau_d = 0.1;
    } else {
            //pos_spherical = pos_to_spherical(p.position[0], p.position[1], p.position[2]);
            pos_scaled = grid_scaling(p.pos_spherical); //position of particle in optical depth grid
            if (pos_scaled[0] == -1.0) {
                //if particle sits outside the grid the optical depth is set to 1e-3
                p.tau_d = 0.001;
            } else {
                //if particle sits inside the grid the optical depth of it is 
                //calculated via a trilinear spline interpolation
                s_phi_values = phi_spline_result(s_phi , radii_v, thetas_v, phis_v, pos_scaled[2]);
                s_theta =  splines_theta ( s_phi_values, radii_v, thetas_v);
                s_theta_values = theta_spline_result( s_theta, radii_v, thetas_v, phis_v , pos_scaled[1]);
                p.tau_d =  tau_p (s_theta_values, radii_v, thetas_v, phis_v, pos_scaled[0]);
                
                if (p.tau_d <0.001) {
                    p.tau_d = 0.001;}
            }
            s_phi_values.clear();
            s_phi_values.shrink_to_fit();
            s_theta.clear();
            s_theta.shrink_to_fit();
            s_theta_values.clear();
            s_theta_values.shrink_to_fit();
        }
  }
  #pragma omp barrier
  for( Particle& p : particles) {
        //update all relevant particle components
        //function RK solver calls the numerical solver to obtain the new particle
        //positions, velocities, size (and mass - consequence of the size change)
        updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
        p.velocity[0], p.velocity[1], p.velocity[2], p.size, p.tau_d, p.temp_d}, current_t, t_next, p.h_updated);
        p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
        p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
        p.size = updated_vector[6];
        p.opac_abs = opac.stellar_abs(p.size);
        p.opac_scat = opac.stellar_scat(p.size);
        p.opac_planck = p.opac_abs + p.opac_scat;
        p.gsca = opac.stellar_gsc(p.size);
        p.mass = dust_mass(p.size); 
        p.h_updated = updated_vector[9];    
        p.temp_d = brent( p.size, p.position[0], 
                                p.position[1], p.position[2], p.tau_d); 
        
        p.pos_spherical = pos_to_spherical(p.position[0], p.position[1], p.position[2]);
        //cout << p.id << endl;
    }
    rm_particles(particles); //removes particles that are too small
    
    if (abs(current_t-plot_time) < 1.0e-8) { 
      long int total = particles.size();
      cout << "Obtaining light curve..." << endl;
      light_curve(particles, current_t);
      if (s_dist == 1) {
      for( Particle& p : particles) {
        output.write((char*) &current_t, sizeof(double));
        output.write((char*) &p.id, sizeof(long int));
        output.write((char*) &p.position[0], sizeof(double));
        output.write((char*) &p.position[1], sizeof(double));
        output.write((char*) &p.position[2], sizeof(double));
        output.write((char*) &p.velocity[0], sizeof(double));
        output.write((char*) &p.velocity[1], sizeof(double));
        output.write((char*) &p.velocity[2], sizeof(double));
        output.write((char*) &p.n_mini, sizeof(double));
        output.write((char*) &p.size, sizeof(double));
        output.write((char*) &p.mass, sizeof(double));
        output.write((char*) &p.tau_d, sizeof(double));
        output.write((char*) &p.temp_d, sizeof(double));
        output.write((char*) &p.opac_planck, sizeof(double));
        output.write((char*) &p.opac_abs, sizeof(double));
        output.write((char*) &p.opac_scat, sizeof(double));
      }
      }
      current_particles = total_particles;
      total_particles = total_particles + nparticles;
      //add particles every 100th of an orbit
      if ((current_t > 0.0) && (abs(current_t-end_t)>1.0e-8) ){
      add_particles(particles, current_particles, total_particles, current_t);
      }
      plot_time = plot_time + major_timestep;
     }

     if (abs(t_next-end_t) < 1.0e-8) {
      cout << "End of simulation.. " << endl;
      vector < dust > dust_grains_out;
      counter = 0;
      for( Particle& p : particles) {
        dust_grains_out.push_back(dust());
        dust_grains_out[counter].timestamp = current_t;
        dust_grains_out[counter].id = p.id;
        dust_grains_out[counter].x_dust = p.position[0];
        dust_grains_out[counter].y_dust = p.position[1];
        dust_grains_out[counter].z_dust = p.position[2];
        dust_grains_out[counter].vx_dust = p.velocity[0];
        dust_grains_out[counter].vy_dust = p.velocity[1];
        dust_grains_out[counter].vz_dust = p.velocity[2];
        dust_grains_out[counter].s_dust = p.size;
        dust_grains_out[counter].nmini = p.n_mini;
        dust_grains_out[counter].h_dust = p.h_updated;
        dust_grains_out[counter].m_dust = p.mass;
        dust_grains_out[counter].temp_dust = p.temp_d;
        dust_grains_out[counter].tau_dust = p.tau_d;
        dust_grains_out[counter].kappa_dust_abs = p.opac_abs;
        dust_grains_out[counter].kappa_dust_scat = p.opac_scat;
        dust_grains_out[counter].kappa_planck = p.opac_planck;
        output_final.write((char *) &dust_grains_out[counter], sizeof(dust));
        counter +=1;
      }

      
      
        memset(extinction, 0.0, sizeof(extinction));
        memset(optical_depth, 0.0, sizeof(optical_depth));
        calculation_ext(particles, extinction, t_global_min);
        optical_depth_calc(extinction, optical_depth);

        for (unsigned long int i = 0; i <t_cells; i++){
            for (unsigned long int j =0; j < p_cells; j++ ) {
                output_tau.write((char*) &thetas_v[i], sizeof(double));
                output_tau.write((char*) &phis_v[j], sizeof(double));
                output_tau.write((char*) &optical_depth[r_cells][i][j], sizeof(double));
            }
        }
      
      }
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<seconds>(stop - start);
     cout << "Iteration took " << duration.count() << " seconds." << endl;
     cout << '\n' << endl;
    current_t = t_next;
    
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
    #pragma omp parallel for
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = cubicInterpolate(p[i], s2, y);
     }
     #pragma omp barrier
    return cubicInterpolate(arr, s1, x);
}

double tricubicInterpolate (vector <vector< vector <double> > > p,
                     vector <double> s1, vector <double> s2, vector <double> s3, 
                     double x, double y, double z) {

    vector <double> arr(s1.size());
    #pragma omp parallel for
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = bicubicInterpolate(p[i], s2, s3, y, z);
     }
    #pragma omp barrier
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



    // if ((current_t >= end_t)) {
        
    //     //memset(extinction, 0.0, sizeof(extinction));
    //     //memset(optical_depth, 0.0, sizeof(optical_depth));
    //     extinction = calculation_ext(particles, t_global_min);
    //     optical_depth = optical_depth_calc(extinction);
       
    //     // optical_depth = optical_depth_calc(extinction);
    //     tau.clear();
    //     tau = tau_to_vector(optical_depth);
    //     s_phi.clear();
    //     //double integrated_tau = 0.0;
    //     s_phi = splines_phi(tau, radii_v, thetas_v, phis_v);

    //     for (unsigned long int i = 0; i <t_cells; i++){
    //         for (unsigned long int j =0; j < p_cells; j++ ) {
    //             output_tau.write((char*) &thetas_v[i], sizeof(double));
    //             output_tau.write((char*) &phis_v[j], sizeof(double));
    //             output_tau.write((char*) &optical_depth[r_cells][i][j], sizeof(double));
    //         }
    //     }
        

    //     for (Particle& p : particles) {
    //         vector <double> pos_spherical, pos_scaled;
    //         pos_spherical = pos_to_spherical(p.position[0], p.position[1], p.position[2]);
    //         pos_scaled = grid_scaling(pos_spherical); //position of particle in optical depth grid
    //         if (pos_scaled[0] == -1.0) {
    //             //if particle sits outside the grid the optical depth is set to 1e-3
    //             p.tau_d = 0.001;
    //         } else {
    //             //if particle sits inside the grid the optical depth of it is 
    //             //calculated via a trilinear spline interpolation
    //             s_phi_values = phi_spline_result(s_phi , radii_v, thetas_v, phis_v, pos_scaled[2]);
    //             s_theta =  splines_theta ( s_phi_values, radii_v, thetas_v);
    //             s_theta_values = theta_spline_result( s_theta, radii_v, thetas_v, phis_v , pos_scaled[1]);
    //             p.tau_d =  tau_p (s_theta_values, radii_v, thetas_v, phis_v, pos_scaled[0]);
                
                
    //             if (p.tau_d < 0.001) {
    //                 p.tau_d = 0.001;}
    //         }
    //         s_phi_values.clear();
    //         s_theta.clear();
    //         s_theta_values.clear();

    //         output_tau_test.write((char*) &p.id, sizeof(long int));
    //         output_tau_test.write((char*) &p.position[0], sizeof(double));
    //         output_tau_test.write((char*) &p.position[1], sizeof(double));
    //         output_tau_test.write((char*) &p.position[2], sizeof(double));
    //         output_tau_test.write((char*) &p.size, sizeof(double));
    //         output_tau_test.write((char*) &p.mass, sizeof(double));
    //         output_tau_test.write((char*) &p.tau_d, sizeof(double));
    //         output_tau_test.write((char*) &p.temp_d, sizeof(double));
    //         output_tau_test.write((char*) &p.opac_planck, sizeof(double));
    //         output_tau_test.write((char*) &p.opac_abs, sizeof(double));
    //         output_tau_test.write((char*) &p.opac_scat, sizeof(double));
    //     }
        
        
        
    // }