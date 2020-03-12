// Defining Variables

#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>

using namespace std;
int i,j,k;

//constants
#define G 6.67408e-11
#define Msun 1.9885e+30
#define Rsun 6.9551e+8
#define AU_to_m 1.496e+11
#define m_to_AU 6.68459e-12

//calculating semi-major axis
double Period_days = 0.8; //period of dust grain in days
double semi;
double Mstar_kg = 0.8*Msun; //mass of star in kg
double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
double semimajor(double period) { //function to evaluate semi-major axis value, period in days, a in AU

   semi = pow(7.496e-6 * Mstar_sun * pow(Period_days, 2.0), (1.0/3.0)); //in AU
   return semi*AU_to_m; //output in meters
}

double a = semimajor(Period_days);

//asign starting indicies
int NR; //number of grid cells
int NP;
int NT;
//int itot = NR + 1; //number of a-mesh indices

//grid boundaries
double rmin;
double rmax;
double Rmin;
double Rmax;

double Pmin;
double Pmax;

double Tmin;
double Tmax;

//grid coordinates (using vectors as arrays) - R, theta, phi
double* Ra; // a-mesh defines edges of grid cells
double* Rb; // b-mesh defines centers of grid cells
double* dRa;
double* dRb;
double Ra_new, Rb_new, dRa_new, dRb_new;

double* Pa; // a-mesh defines edges of grid cells
double* Pb; // b-mesh defines centers of grid cells
double* dPa;
double* dPb;
double Pa_new, Pb_new, dPa_new, dPb_new;

double* Ta; // a-mesh defines edges of grid cells
double* Tb; // b-mesh defines centers of grid cells
double* dTa;
double* dTb;
double Ta_new, Tb_new, dTa_new, dTb_new;

//gaussian grid
//R grid
double* R_x;
double* R_g;
double* R_inv;
double* DR;
double R_A, R_B, R_C; //scale factor of gaussian - sets max DR
double R_mu, R_sd;
double R_suminv, sumDR;

//PHI grid
double* P_x;
double* P_g;
double* P_inv;
double* DP;
double P_A, P_B, P_C; //scale factor of gaussian - sets max DR
double P_mu, P_sd;
double P_suminv, sumDP;

//THETA grid
double* T_x;
double* T_g;
double* T_inv;
double* DT;
double T_A, T_B, T_C; //scale factor of gaussian - sets max DR
double T_mu, T_sd;
double T_suminv, sumDT;

//density variables
double density = 1.83e-12; //kg/m3 (av density)
double density_bulk = 3000.; //kg/m3 (bulk density)
double T_mean;
double T_stde;
double P_mean;
double P_stde;
double R_mean;
double R_stde;

double mean;
double stde;

//initializing density array
double*** d;
double d_new_R, d_new_P, d_new_T;
double*** total_mass;
double TMASS;
double tmass;

//opacity variables
double*** kappa; //dimensionless opacity

//optical depth variables
double*** t_num; //optical depth (dimensionless anyway)
double*** t_ana;
double*** t;

//counter thing
double* T_vec;
double* P_vec;
double* R_vec;
double*** den;
double mass; //mass of particles
int noparticles;
int nop;
