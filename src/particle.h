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
    vector <double> pos_p; //firt projection for lc calculation
    vector <double> pos_dp; //second projection for lc calculation
    double size; //size of the particle
    double n_mini; //number of grains inside superparticle
    double opac_abs;
    double opac_scat; 
    double opac_planck;
    double gsca;
    double h_updated; //current optimal time step for particle
    double mass; //mass of particle
    double tau_d; //optical depth
    double temp_d; //particle temperature
    double f_scat;
    double err;

};

void add_particles(vector <Particle> &particles,long int &current_particles, long int &total_particles, double total_t);



void solve_particles(double total_t, double end_t, vector <Particle> &particles, \
                     long int total_particles,long int current_particles);

void rm_particles(vector <Particle> &particles);

void calculation_ext(vector <Particle>& particles, double (&ext)[r_cells][t_cells][p_cells],  double delta_t);

//light curve

void light_curve(vector<Particle>& particles, double current_t);
void  extinction_lc( vector <Particle>& particles, vector <vector <double>> &patches, 
         vector<double> &h_grid, vector<double> &v_grid,
         vector<vector <double>> &taus, int h_cells, int v_cells,
         double z_max, double z_min, double y_max, double y_min, double r_star_a);


extern vector <Particle> particles;
extern bool tau_constant;
