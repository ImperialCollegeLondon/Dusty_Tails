#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "particle.h"
#include "functions.h"
#include "opacities.h"
#include <omp.h>



using namespace std;
vector <double> rad_vector(3);

vector <double> centrifugal(vector <double> V){
  //function to calculate centrifugal force
  vector <double> omxr, finalx(3);
  //cout << ang_vel << endl;
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
  // and to evaluate poynting roberston drag
	vector <double> v_drag(3), pr_vector(3);
  double beta, kappa, constant;
  double d_prod;
  vector <double> rad_vector(3);
  //cout << V[6] << endl;
  kappa =opac.stellar_abs(V[6]) + opac.stellar_scat(V[6]);
  //cout << kappa << endl;
  beta = beta_fn(kappa, V[7], V[6]);
  //cout << "previous kappa " << (3.0/4.0) / (rho_d * 1.e-4) << endl ;
  //cout << "current kappa " << kappa << endl;

  constant = (beta*G_dim)/(pow(scalar(V[0]- star_pos[0], V[1], V[2]), 3.0));
  
  d_prod = dot_product({V[3], V[4], V[5]}, sunit_vector(V));
  v_drag = drag_vel(V);

  for (unsigned int i = 0; i < 3; i++) {
      rad_vector[i] = constant* (1-((d_prod)/c_dim)) *(V[i]-star_pos[i]) - constant*v_drag[i]/c_dim;
 }
  return rad_vector;

}
