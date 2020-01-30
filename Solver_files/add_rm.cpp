#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"


using namespace std;


ofstream file("output5.txt", std::ios_base::out | std::ios_base::app);

vector <Particle> add_particles(vector <Particle> particles, long int current,
                   long int total, double time){

        double v_esc;

        v_esc = (pow((2.0*G *0.1*Mearth)/(0.5*Rearth), 0.5)) * (T/a);
        //escape velocity m/s

        for ( unsigned long int i = current; i < total; i++){

                 Particle grain;
                 grain.id = i+1;
                 double theta = fRand(0.0, PI);

                 double phi = fRand(0.0, 2.0*PI);


                 grain.position = {r_start*sin(theta)*cos(phi) + planet_x, \
                                   r_start*sin(theta)*sin(phi), \
                                   r_start*r_planet_dim*cos(theta)};

                 grain.velocity = {v_esc*sin(theta)*cos(phi), \
                                   v_esc*sin(theta)*sin(phi), \
                                   v_esc*cos(theta)};
                 grain.p_size = size;
                 grain.p_tau = tau;
                 grain.p_density = rho_d;
                 grain.h_updated = 0.001;

                 particles.push_back(grain);

                 file << time << ",";
                 file << grain.id << ",";
                 file << grain.position[0] << ",";
                 file << grain.position[1] << ",";
                 file << grain.position[2] << "\n";
               }


              return particles;


}

void solve_particles(double total_t, double end_t, vector <Particle> particles, \
                     long int total_particles, double t_common, double big_step, \

                    long int current_particles){

  double plot_time = 0.5;

  while (total_t < end_t) {

    for( Particle& p : particles) {
      vector <double> updated_vector;
      updated_vector.clear();
      updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
      p.velocity[0], p.velocity[1], p.velocity[2]}, total_t, t_common, p.h_updated);

      p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
      p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};

      p.h_updated = updated_vector[6];

      file << total_t + big_step << ",";
      file << p.id << ",";
      file << p.position[0] << ",";
      file << p.position[1] << ",";
      file << p.position[2] << "\n";


    }



    if (total_t > plot_time) {
      current_particles = total_particles;
      total_particles = total_particles + 10;
      particles = add_particles(particles, current_particles, total_particles, total_t);
      plot_time = plot_time + 0.5;
    }

    total_t = total_t + big_step;
    t_common = t_common + big_step;


  }
}
