#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "particle.h"
#include "functions.h"


using namespace std;
vector <double> rad_vector(3);

vector <double> centrifugal(vector <double> V){
  //function to calculate centrifugal force
  vector <double> omxr, finalx(3);
  omxr = cross_product(0.0, 0.0, ang_vel, V[0], V[1], V[2]);
  finalx = cross_product(0.0, 0.0, ang_vel, omxr[0], omxr[1], omxr[2]);

  return finalx;
}

vector <double> coriolis(vector <double> V){
  //function to calculate coriolis force
  vector <double> omxv,  rvxr, finalx;
  omxv = cross_product(0.0, 0.0, ang_vel, V[3], V[4], V[5]);
  return omxv;

}

vector <double> rad_pressure(vector <double> V){
  //function to calculate radiation pressure force, accounting for red shift
  double beta, k, constant;
  double d_prod;
  vector <double> rad_vector(3);

  k = opacity(V[6], V[0], V[1], V[2]);
  beta = beta_fn(k, V[7]);


  constant = (beta*G_dim)/(pow(scalar(V[0]- star_pos[0], V[1], V[2]), 3.0));
  if (isnan(k)) {
    cout << "oops opacity is nan " << endl;
  }

  d_prod = dot_product({V[3], V[4], V[5]}, sunit_vector(V));
  

  for (unsigned int i = 0; i < 3; i++) {

      rad_vector[i] = constant* (1-((d_prod)/c_dim)) *(V[i]-star_pos[i]);
 }

  return rad_vector;

}

vector <double> pr_drag(vector <double> V){
  //function to evaluate poynting roberston drag
	vector <double> v_drag(3), pr_vector(3);

	double constant, beta, k, Lum;

  k = opacity(V[6], V[0], V[1], V[2]);
  beta = beta_fn(k, V[7]);

  constant = (beta*G_dim)/(pow(scalar(V[0]-star_pos[0], V[1], V[2]), 3.0)*c_dim);

  v_drag = drag_vel(V);

  for (unsigned int i=0; i < 3; i++) {
      pr_vector[i] = constant*v_drag[i];
  }
  return pr_vector;
}
