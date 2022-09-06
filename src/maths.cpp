#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "particle.h"
#include "functions.h"
#include <omp.h>


using namespace std;


double scalar(double x, double y, double z){
        double s =0.0;
        s = pow( pow(x, 2.) + pow(y, 2.) + pow(z, 2.), 0.5 );
        return s;
}

vector <double> cross_product(double m1, double m2, double m3, double n1, double n2, double n3){
  //function to calculate the cross product between two vectors
  vector <double> cp_vector(3, 0.0);

	double i_new =0.0, j_new=0.0, k_new=0.0;

	i_new = m2*n3 - m3*n2;
	j_new = m3*n1 - m1*n3;
	k_new = m1*n2 - m2*n1;

  cp_vector = {i_new, j_new, k_new};

	return cp_vector;


}

double dot_product(vector <double> n,  vector <double> m){
     //function to calculate the dot product between two vectors
     double dp = (n[0]*m[0]) + (n[1]*m[1]) + (n[2]*m[2]);

		 return dp;

}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / (double)RAND_MAX;
    return fMin + f * (fMax - fMin);
}

vector <double> pos_to_spherical(double x, double y, double z){
  double x_new, radius, theta, phi;
  vector <double> s_pos(3, 0.0);
  x_new = x - planet_x + 1.0;

  radius = pow( pow(x_new, 2.) + pow(y, 2.) + pow(z, 2.), 1./2.);
  theta = acos(z / radius);
  phi = atan(y / x_new);

  s_pos = {radius, theta, phi};

  return s_pos;

}

vector <double> vel_to_spherical(double x, double y, double z){
  double v_radius, v_theta, v_phi;
  vector <double> v_pos(3, 0.0);

  v_radius = pow( pow(x, 2.) + pow(y, 2.) + pow(z, 2.), 1./2.);
  v_theta = acos(z / v_radius);
  v_phi = atan(y / x);

  v_pos = {v_radius, v_theta, v_phi};

  return v_pos;

}
