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
double NR; //number of grid cells
//int itot = NR + 1; //number of a-mesh indices

//grid boundaries
double rmin = 0.;
double rmax = 2.*a;
double Rmin=rmin/a;
double Rmax=rmax/a;

//grid coordinates (using vectors as arrays)
vector <double> Ra; // a-mesh defines edges of grid cells
vector <double> Rb; // b-mesh defines centers of grid cells
vector <double> dRa;
vector <double> dRb;
double Ra_new, Rb_new, dRa_new, dRb_new;

//density variables
double density = 5000.; //kg/m3 (av density)
vector <double> d; //dimensionless density
double mean;
double sd;
vector <double> gauss;
double gauss_new;

//opacity variables
vector <double> k; //dimensionless opacity
double k_new; //used to fill opacity

//optical depth variables
vector <double> t; //optical depth (dimensionless anyway)
double t_new;
