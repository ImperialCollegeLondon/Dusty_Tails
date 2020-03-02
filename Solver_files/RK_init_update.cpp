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

vector <Particle> particles;

int main() {

  long int total_particles = 50;
  double t_common = 0.01; //arbitrary for now
  double big_step = 0.01;
  double end_t = 2.0;
  double total_t = 0.0;

  //srand(2);

  long int current_particles = 0;

  particles = add_particles(particles, current_particles, total_particles, 0.0);

  particles = rm_particles(particles);

  solve_particles(total_t, end_t, particles, total_particles, t_common, big_step,
                       current_particles);


}
