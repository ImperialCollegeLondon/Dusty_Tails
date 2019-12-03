//Header file with constants


#define PI 3.14159

#define G 6.67408e-11
#define Msun 1.9885e+30
#define Rsun 6.9551e+8
#define AU_to_m 1.496e+11
#define m_to_AU 6.68459e-12
#define Mearth 5.972e+24
#define c 2.9979245e+8


const double Mstar_kg = 0.8*Msun; //mass of star in kg
const double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
const double Period_days = 0.8; //period of dust grain in days
const double T = 0.8*24.0*60.0*60.0; //period of dust grain in seconds

const double Temp = 4400.0; //stars temperature
const double S = 0.9; //safety factor
const double tol = 1e-6; //error tolerance

const double beta = 0.0; //ratio of radiation pressure to gravity

const double tau = 0.1; //optical depth

//constants in cgs units for beta calculation

#define sigma 5.6705119e-5 //stefan boltzmann
#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define rho_d 3.0 //density of dust particle
#define Msun_cgs 1.9885e+33 // solar mass
#define s 1.0e-4 //dust particle size
