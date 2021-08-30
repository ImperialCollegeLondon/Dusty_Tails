#include <iostream>
#include <vector>
#include "spline.h"
#include <array>

using namespace std;
const int r_cells = 200;
const double r_cells_d = 200.;
const int t_cells = 50;
const double t_cells_d = 50.;
const int p_cells = 250;
const double p_cells_d = 250.;

const double d_r_min = 0.9;
const double d_r_max = 1.1;
const double d_t_min = 1.55;
const double d_t_max = 1.60;
const double d_p_min = -0.25;
const double d_p_max = 0.0;

const double d_dr = (d_r_max - d_r_min)/ r_cells_d;
const double d_dtheta = (d_t_max - d_t_min ) / t_cells_d;
const double d_dphi = ( d_p_max - d_p_min) / p_cells_d;

void error(double analytic [r_cells][t_cells][p_cells], double numerical [r_cells][t_cells][p_cells], double errors [r_cells][t_cells][p_cells]);

extern double r_a [r_cells + 1];
extern double r_b [r_cells];
extern double theta_a [t_cells + 1];
extern double theta_b [t_cells];
extern double phi_a [p_cells + 1];
extern double phi_b [p_cells];

extern vector <double> radii_v, thetas_v, phis_v;

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
extern vector <Particle> particles;

void add_particles(vector <Particle> &particles, long int current_particles,
                   long int total_particles, double total_t);


void solve_particles(double t_global, double end_t, vector <Particle> &particles, \
                     long int total_particles,long int current_particles);

void rm_particles(vector <Particle> &particles);
vector <vector < vector <double> > > calculation_ext();
vector <vector < vector <double> > > optical_depth_calc(vector <vector < vector <double> > > ext);

extern vector <Particle> particles;

vector <double> r_grid_to_vector(double r[r_cells+1]);
vector <double> t_grid_to_vector(double t[t_cells+1]);
vector <double> p_grid_to_vector(double p[p_cells+1]);

void test_pos(double *r_test, double *theta_test, double *phi_test, \
                double d_r_min, double d_r_max, double d_t_min, double d_t_max,
                double d_p_min, double d_p_max, double d_dr, double d_dtheta, double d_dphi);

void test_dist_ext(double *r_test, double *theta_test, double *phi_test, double function[r_cells][t_cells][p_cells]);
vector < vector < vector <double> > >  tau_to_vector(vector < vector < vector <double> > > tau);