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

int main() {
  double init_vel = pow((G_dim*(m_planet)/(0.1*r_h)), 0.5); //inertial frame

  vector < Particle > particles;

  vector < vector <double> > initial_pos;
  vector < vector <double> > initial_vel;

  vector <double> pos_1;
  vector <double> pos_2;
  vector <double> vel_1;
  vector <double> vel_2;

  double v_esc;

  v_esc = (pow((2.0*G *0.1*Mearth)/(0.5*Rearth), 0.5)) * (T/a); //escape velocity m/s
  pos_1 = {planet_x, 0.0, 0.0};
  pos_2 = {planet_x, 0.0, 0.0};


  vel_1 = {v_esc, 0.0, 0.0};
  vel_2 = {-v_esc, 0.0, 0.0};

  initial_pos.push_back(pos_1);
  initial_pos.push_back(pos_2);

  initial_vel.push_back(vel_1);
  initial_vel.push_back(vel_2);

  long int no_particles = 5;
  double t_common = 0.1; //arbitrary for now
  double big_step = 0.1;
  double end_t = 10.0;
  double total_t = 0.0;

  long int current_particles = 0;


  ofstream file("data1.txt");

  for ( unsigned long int i = current_particles; i < no_particles; i++){
      Particle grain;
      grain.id = i;
      grain.position = initial_pos[i];
      grain.velocity = initial_vel[i];
      grain.p_size = size;
      grain.p_tau = tau;
      grain.p_density = rho_d;
      grain.h_updated = 0.001;

      particles.push_back(grain);
      file << total_t << ",";
      file << grain.id << ",";
      file << grain.position[0] << ",";
      file << grain.position[1] << ",";
      file << grain.position[2] << "\n";

      current_particles = current_particles + 1;
    }



  while (total_t < end_t) {


    for( Particle& p : particles) {
      vector <double> updated_vector;
      updated_vector.clear();
      updated_vector = RK_solver({p.position[0], p.position[1], p.position[2], \
      p.velocity[0], p.velocity[1], p.velocity[2]}, total_t, t_common, p.h_updated);
      p.position = {updated_vector[0],updated_vector[1], updated_vector[2]};
      p.velocity = {updated_vector[3],updated_vector[4], updated_vector[5]};
      file << total_t + big_step << ",";
      file << p.id << ",";
      file << p.position[0] << ",";
      file << p.position[1] << ",";
      file << p.position[2] << "\n";

    }


    total_t = total_t + big_step;
    t_common = t_common + big_step;

    if (total_t

  }

}
