//Header file with constants


#define PI 3.14159

#define G 6.67408e-11
#define Msun 1.9885e+30
#define Rsun 6.9551e+8
#define AU_to_m 1.496e+11
#define m_to_AU 6.68459e-12

const double Mstar_kg = 0.8*Msun; //mass of star in kg
const double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
const double Period_days = 0.8; //period of dust grain in days
const double T = 0.8*24.0*60.0*60.0; //period of dust grain in seconds

const double S = 0.9; //safety factor
const double tol = 1e-7; //error tolerance