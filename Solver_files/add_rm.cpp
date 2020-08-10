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
mt19937 generator (seed);
uniform_real_distribution<double> uniform_phi(0.625, 0.875);
uniform_real_distribution<double> uniform_theta(0.2, 0.8);

//uniform_real_distribution<double> uniform_phi(0.0, 1.0);
//uniform_real_distribution<double> uniform_theta(0.0, 1.0);

ofstream ofile("kic_k222_035.bin", ios::out | ios::binary);

void add_particles(vector <Particle> &particles, long int current,
                   long int total, double time){

        double v_esc;
        Particle grain;

        v_esc = (pow((2.0*G *0.05*Mearth)/(0.38*Rearth), 0.5)) * (T/a);
        //escape velocity m/s

        for ( unsigned long int i = current; i < total; i++){


                 grain.id = i+1;
                 double phi = 2.0*PI * uniform_phi(generator);
                 //double phi = 2.0*PI * 0.7;

                 double theta = acos(1- 2.0* uniform_theta(generator));
                 //double theta = acos(1- 2.0* 0.55);


                 grain.position = {r_start*sin(theta)*cos(phi) + planet_x, \
                                   r_start*sin(theta)*sin(phi), \
                                   r_start*cos(theta)};

                 grain.velocity = {v_esc*sin(theta)*cos(phi), \
                                   v_esc*sin(theta)*sin(phi), \
                                   v_esc*cos(theta)};

                 grain.p_size = 0.35e-4;
                 grain.p_tau = tau;
                 grain.p_density = rho_d;
                 grain.h_updated = 0.001;
                 grain.p_mass = dust_mass(grain.p_size);

                 particles.push_back(grain);

               }

}

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

        ofile.write((char*) &total_t, sizeof(double));

        ofile.write((char*) &p.id, sizeof(long int));
        cout << "ID " << p.id << endl;

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

    if (total_t > plot_time) {
      current_particles = total_particles;
      total_particles = total_particles + 100;
      add_particles(particles, current_particles, total_particles, total_t);
      plot_time = plot_time + 0.01;
  }


    total_t = total_t + big_step;
    t_common = t_common + big_step;


  }
}
