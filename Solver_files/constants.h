//Header file with constants
#include <iostream>
#include <cmath>
#include <vector>


using namespace std;
#define PI 3.14159
#define wien 0.289 //Wiens constant in cgs
#define G 6.67408e-11
#define Msun 1.9885e+30
#define Rsun 6.9551e+8
#define AU_to_m 1.496e+11
#define m_to_AU 6.68459e-12
#define Mearth 5.972e+24
#define c 2.9979245e+8
#define Rearth 6.3781e+6

const double Mstar_kg = 0.67*Msun; //mass of star in kg
const double Mstar_sun = 0.67; //mass of star in terms of mass of the sun
const double Period_days = 15.68/24.0; //period of planet in days
const double T = 15.68*60.0*60.0; //period of planet in seconds
const double Rstar = 0.67;


const double a = pow((G*Mstar_kg* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0);
const double m_planet = (0.05*Mearth)/Mstar_kg;
const double r_planet = 0.38*Rearth;
const double r_start = (2.*r_planet)/a;
const double r_planet_dim = r_planet/a;
const double G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant
const double c_dim = c * (T / a);
const double r_h = pow(m_planet/3.0, 1.0/3.0); //hill radius

extern double star_x;
extern double planet_x;

extern vector <double> star_pos;

extern vector <double> planet_pos;


//const double h0 = 0.001; //initial time step

const double Temp = 4550.0; //stars temperature
const double S = 0.9; //safety factor
const double tol = 1e-8; //error tolerance
const double A = 77365.; //clausius claperyon relation
const double B = 39.3; //clausius claperyon relation
const double alpha = 0.1;
const double tau = 0.1; //optical depth
const double mu = 101.961; //molecular weight

//constants in cgs units for beta calculation
#define amu 1.661e-24
#define kb 1.381e-16
#define sigma 5.6705119e-5 //stefan boltzmann
#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define rho_d 4.0 //density of dust particle
#define Msun_cgs 1.9885e+33 // solar mass
 //dust particle size
