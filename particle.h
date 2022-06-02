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
    double size; //size of the particle
    double opac_abs;
    double opac_scat; 
    double opac_planck;
    double h_updated; //current optimal time step for particle
    double mass; //mass of particle
    double tau_d; //optical depth
    double temp_d; //particle temperature

};

void add_particles(vector <Particle> &particles,long int &current_particles, long int &total_particles, double total_t);



void solve_particles(double total_t, double end_t, vector <Particle> &particles, \
                     long int total_particles,long int current_particles);

void rm_particles(vector <Particle> &particles);

double*** calculation_ext(vector <Particle>& particles,  double delta_t);

extern vector <Particle> particles;
extern bool tau_constant;
