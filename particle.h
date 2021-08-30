#include <iostream>
#include <vector>
#include "spline.h"

using namespace std;

class Particle {
   public:
    long int id; //unique id of the particle
    vector <double> position; //current position of the particle
    vector <double> velocity; //current velocity of the particle
    vector <double> pos_spherical; // current position in spherical coordinates
    vector <double> v_spherical; //curretn velocity in spherical coordinates
    double p_size; //size of the particle
    double p_density; //particle density
    double h_updated; //current optimal time step for particle
    double p_mass; //mass of particle
    double p_tau; //optical depth

};

void add_particles(vector <Particle> &particles, long int current_particles,
                   long int total_particles, double total_t);



void solve_particles(double total_t, double end_t, vector <Particle> &particles, \
                     long int total_particles,long int current_particles);

void rm_particles(vector <Particle> &particles);

void calculation_ext(vector <Particle>& particles, double ext [r_cells][t_cells][p_cells]);

extern vector <Particle> particles;
