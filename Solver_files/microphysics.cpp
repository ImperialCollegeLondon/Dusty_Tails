#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"


using namespace std;

double omega(double mplanet, double mstar){
    return pow((G_dim *(m_planet + mstar)), 0.5);
}

double beta_fn(double k, double L_star, double M_star){
  //mass of star in terms of mass of sun, and everything in cgs units
  //function to evaluate beta (radiation accel/ gravity accel)
  double beta =0.0;
	beta = (k*L_star*exp(-tau)) / (4.0*PI*c_cgs*G_cgs*M_star*Msun_cgs);
	return beta;

}

double opacity(double s, double x, double y, double z){
    //function to calculate opacity in cgs units
    //might need to add density and size as parameters later
    //cout << "Opacity " << (3.0/4.0)* (qfactor(s, x, y, z)/ (rho_d * s)) << endl;
    return (3.0/4.0)* (qfactor(s, x, y, z)/ (rho_d * s)) ;
}

double qfactor(double s, double x, double y, double z){
  //cgs units
  return (s*temp_dust(luminosity(Rstar),s, x, y, z)) / wien;
}

double clausius_clap(double s, double x, double y, double z){
  double cc =0.0;

  cc = exp((-A/temp_dust(luminosity(Rstar), s, x, y, z)) + B);

  return cc;

}

double luminosity(double R_star){
  //function to evaluate stars luminosity with stefan boltzmann law
  //cgs units
	return sigma*4.0*PI* pow(R_star*Rsun_cgs, 2.0) * pow(Temp, 4.0);

}

double radial_vel(vector <double> vel, vector <double> s_vector){
  //function to obtain radial velocity term
	return dot_product(vel, s_vector);

}

vector <double> drag_vel(vector <double> V){
   //function to obtain velocity relative to radiation field
	 vector <double> omegar, vrad(3);

     double ang_vel = omega(m_planet, 1.0);

	 omegar = cross_product(0.0, 0.0, ang_vel, V[0], V[1], V[2]);

     for (unsigned int i=0; i < 3; i++){
         vrad[i] = V[i+3] + omegar[i];
    }
	 return vrad;
}

vector <double> sunit_vector(vector <double> V){
  vector <double> s_unit(3);

  for (unsigned int i = 0; i < 3; i++) {
  //function to obtain unit vector of radiation field
    s_unit[i] = (V[i] - star_pos[i]) / scalar(V[0] - star_pos[0],V[1],V[2]);

   }

	return s_unit;
}

double temp_dust(double lum, double s, double x, double y, double z){
  double Tdust, dl;
  dl = scalar((x-star_pos[0]), y, z)*a*pow(10., 2.0);

  Tdust = pow((lum*wien)/(16.0*PI*pow(dl, 2.0)*sigma*s), 1./5.);

  return Tdust;
}

double dust_mass(double s){
    double md;
    md = rho_d * (4.0/3.0) * PI * pow(s, 3.0);
    return md;
}
