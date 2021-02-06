//File that has the main function

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include <stdlib.h>
#include <time.h>

using namespace std;

vector <Particle> particles; //initiate vector of "Particle" (Object defined in particle.h)

int main() {

  long int total_particles = 1000; //initial number of particles to start simulation with
  double t_common = 0.01;
  double big_step = 0.01; //big time step (in terms of planetary orbits)
  double end_t = 1.0; // end time of simulation
  double total_t = 0.0; // total time that has passed, so 0 in the beginning

  long int current_particles = 0; // number of current particles in simulation

  add_particles(particles, current_particles, total_particles, 0.0); // call function that adds particles to simulation (in add_rm.cpp file)

  rm_particles(particles); //remove any particles that were generated too far from the planet

  solve_particles(total_t, end_t, particles, total_particles, t_common, big_step,
                       current_particles); //solver function (in add_rm.cpp file)


}
