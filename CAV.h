// Defining Variables

#include <vector>
#include <cmath>

using namespace std;
int i;

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
double NR, NT, NP; //number of grid cells
//int itot = NR + 1; //number of a-mesh indices

//grid boundaries
double rmin = 0.;
double rmax = 2.*a;
double Rmin=rmin/a;
double Rmax=rmax/a;

double Tmin = 0.;
double Tmax = math.pi;

double Pmin = 0.;
double Pmax = (math.pi)/2.;

//grid coordinates (using vectors as arrays) - R, theta, phi
vector <double> Ra; // a-mesh defines edges of grid cells
vector <double> Rb; // b-mesh defines centers of grid cells
vector <double> dRa;
vector <double> dRb;
double Ra_new, Rb_new, dRa_new, dRb_new;

vector <double> Ta;
vector <double> Tb;
vector <double> dTa;
vector <double> dTb;
double Ta_new, Ta_new, dTa_new, dTb_new;

vector <double> Pa;
vector <double> Pb;
vector <double> dPa;
vector <double> dPb;
double Pa_new, Pb_new, dPa_new, dPbs_new;

//gaussian grid
//R grid
vector <double> g_R,DR,x_R,inv_R;
double A_R, B_R, C_R; //scale factor of gaussian - sets max DR
double B_R; //addition factor of gaussian - sets min DR
double mu_R, sd_R;
double suminv_R, inv_new_R, sumDR, nr, g_new_R;

//Theta grid
vector <double> g_T,DT,x_T,inv_T;
double A_T, B_T, C_T; //scale factor of gaussian - sets max DR
double B_T; //addition factor of gaussian - sets min DR
double mu_T, sd_T;
double suminv_T, inv_new_T, sumDT, nt, g_new_T;

//Phi grid
vector <double> g_P,DP,x_P,inv_P;
double A_P, B_P, C_P; //scale factor of gaussian - sets max DR
double B_P; //addition factor of gaussian - sets min DR
double mu_P, sd_P;
double suminv_P, inv_new_P, sumDP, np, g_new_P;

//density variables
double density = 1.83e-12; //kg/m3 (av density)
vector <double> d; //dimensionless density
double mean;
double stde;
vector <double> gauss;
double gauss_new;
double density_bulk = 3000.; //kg/m3 (bulk density)

//opacity variables
vector <double> k; //dimensionless opacity
double k_new; //used to fill opacity

//optical depth variables
vector <double> t; //optical depth (dimensionless anyway)
double t_new;
