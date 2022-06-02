//Header file with constants
#include <iostream>
#include <cmath>
#include <vector>
#include "spline.h"


using namespace std;

#define PI 3.14159
#define wien 0.289 //Wiens constant in cgs
#define G 6.67408e-11 //gravitational constant
#define Msun 1.9885e+30 //solar mass in kg
#define Rsun 6.9551e+8 //solar radius in meters
#define AU_to_m 1.496e+11 //AU in meters
#define m_to_AU 6.68459e-12 //meters in AU
#define Mearth 5.972e+24 //Earth mass in kg
#define clight 2.9979245e+8 //speed of light in SI
#define Rearth 6.3781e+6 //Earth radius in meters

//definition of the grid cell numbers for the optical depth grid
const int r_cells = 96;
const double r_cells_d = 96.;
const int t_cells = 40;
const double t_cells_d = 40.;
const int p_cells =120;
const double p_cells_d = 120.;

//declaration of some external variables
extern double star_x; //x position of the star, w.r.t CoM of system
extern double planet_x; //x position of the planet, w.r.t CoM of system
extern vector <double> star_pos; //vector of the stars position
extern vector <double> planet_pos; //vector of the planets position

//Frame of reference parameters
extern double ang_vel;

//Stellar parameters:
extern  double Mstar_kg, Mstar_sun, Rstar, Temp, lum; 

//Planetary parameters:
extern  double Period_days, T, a, m_planet, r_planet, r_start, r_planet_dim, r_h;

//Dust parameters:
extern  double A, Bp, rho_d, s_0,alpha, mu;
extern  string opac_data, opacity_dir;
extern  string comp, outflow_s, T_int_s, output_file;

//Outflow parameters:
extern double mdot, v_esc, mdot_read;
extern int outflow, tau_type;
extern bool tau_constant;
extern int cont;

//Grid parameters:
extern double d_r_min, d_r_max,  d_t_min, d_t_max, d_p_min, d_p_max;

extern long int current_particles, total_particles; // number of current particles in simulation
extern double dr, dtheta, dphi;
extern double d_dr, d_dtheta, d_dphi;

//Some dimensionless quantitites:
extern double G_dim, c_dim;

//Numerical factors for integrator:
const double S = 0.9; //safety factor
const double tol = 1.0e-7; //error tolerance

//GYR in seconds
const double gyr = pow(10.,9) * 365. * 24. * 60. * 60. ;// 1Gyr in seconds


//constants in cgs units for beta calculation
#define amu 1.661e-24
#define kb 1.381e-16
#define sigma 5.6705119e-5 //stefan boltzmann
#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define Msun_cgs 1.9885e+33 // solar mass
#define Mearth_cgs 5.972e+27 // earth mass grams
