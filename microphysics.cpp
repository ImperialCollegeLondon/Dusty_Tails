#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
//#include "RK_variables.h"
#include "functions.h"

using namespace std;

double omega(double mplanet, double mstar){
    return pow((G_dim *(m_planet + mstar)), 0.5);
}

double beta_fn(double k, double L_star, double M_star){
  //mass of star in terms of mass of sun, and everything in cgs units
  //function to evaluate beta (radiation accel/ gravity accel)
  double beta;
	beta = (k*L_star*exp(-tau)) / (4.0*PI*c_cgs*G_cgs*M_star*Msun_cgs);
	return beta;

}

double opacity(double Q_fn){
    //function to calculate opacity in cgs units
    //might need to add density and size as parameters later
    return (3.0/4.0)* (Q_fn/ (rho_d * size)) ;
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

vector <double> drag_vel(double x, double y, double z, double vx, double vy, double vz){
   //function to obtain velocity relative to radiation field
	 vector <double> omegar, vrad;
	 vrad.clear();

	 omegar = cross_product(0.0, 0.0, 2.0*PI, x, y, z);

	 vrad.push_back(vx + omegar[0]);
	 vrad.push_back(vy + omegar[1]);
	 vrad.push_back(vz + omegar[2]);

	 return vrad;
}

vector <double> s_vector(double x, double y, double z){
  vector <double> s_unit;
  //function to obtain unit vector of radiation field
  s_unit.clear();
	s_unit.push_back((x - star_pos[0]) /scalar(x,y,z));
	s_unit.push_back((y - star_pos[1]) /scalar(x,y,z));
	s_unit.push_back((z - star_pos[2]) /scalar(x,y,z));

	return s_unit;
}
