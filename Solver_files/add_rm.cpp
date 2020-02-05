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

ofstream ofile("output3.bin", ios::out | ios::binary);

vector <Particle> add_particles(vector <Particle> particles, long int current,
                   long int total, double time){

        double v_esc;

        v_esc = (pow((2.0*G *0.05*Mearth)/(0.38*Rearth), 0.5)) * (T/a);
        //escape velocity m/s

        for ( unsigned long int i = current; i < total; i++){

                 Particle grain;
                 grain.id = i+1;
                 double theta = fRand(0.0, PI);

                 double phi = fRand(0.0, 2.0*PI);

                 /*
                 grain.position = {r_start*sin(theta)*cos(phi) + planet_x, \
                                   r_start*sin(theta)*sin(phi), \
                                   r_start*r_planet_dim*cos(theta)};

                 grain.velocity = {v_esc*sin(theta)*cos(phi), \
                                   v_esc*sin(theta)*sin(phi), \
                                   v_esc*cos(theta)}; */
                 grain.position = {planet_x, 0.0, 0.0};
                 grain.velocity = {v_esc, 0.0, 0.0};
                 grain.p_size = dsize;
                 grain.p_tau = tau;
                 grain.p_density = rho_d;
                 grain.h_updated = 0.001;

                 grain.p_temp = temp_threshold(1.0, luminosity(Rstar),\
                                 grain.position[0], grain.position[1],
                                 grain.position[2]);

                 particles.push_back(grain);


                 ofile.write((char*) &time, sizeof(double));
                 ofile.write((char*) &grain.id, sizeof(long int));
                 ofile.write((char*) &grain.position[0], sizeof(double));
                 ofile.write((char*) &grain.position[1], sizeof(double));
                 ofile.write((char*) &grain.position[2], sizeof(double));
                 ofile.write((char*) &grain.p_temp, sizeof(double));

               }


              return particles;


}

vector <Particle> rm_particles(vector <Particle> particles){
  for (unsigned long i = 0; i < particles.size(); i++){
    if (particles[i].p_temp > 2700.0) {
      particles.erase(particles.begin() + i);
      i--;
    }
  for (unsigned long i = 0; i < particles.size(); i++){
    if (fabs(particles[i].position[0]) > 3.0) {
        particles.erase(particles.begin() + i);
        i--;
    }
  }
  for (unsigned long i = 0; i < particles.size(); i++){
    if (fabs(particles[i].position[1]) > 3.0) {
        particles.erase(particles.begin() + i);
        i--;
    }
   }

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

      p.p_temp = temp_threshold(1.0, luminosity(Rstar),\
                      p.position[0], p.position[1],
                      p.position[2]);

      double time_now = total_t + big_step;

      ofile.write((char*) &time_now, sizeof(double));
      ofile.write((char*) &p.id, sizeof(long int));
      ofile.write((char*) &p.position[0], sizeof(double));
      cout << "x " << p.position[0] << endl;
      ofile.write((char*) &p.position[1], sizeof(double));
      cout << "y " << p.position[1] << endl;
      ofile.write((char*) &p.position[2], sizeof(double));
      cout << "z " << p.position[2] << endl;
      ofile.write((char*) &p.p_temp, sizeof(double));



    }

    //particles = rm_particles(particles);
    /*
    if (total_t > plot_time) {
      current_particles = total_particles;
      total_particles = total_particles + 500;
      particles = add_particles(particles, current_particles, total_particles, total_t);
      plot_time = plot_time + 0.5;
  } */


    total_t = total_t + big_step;
    t_common = t_common + big_step;


  }
}
