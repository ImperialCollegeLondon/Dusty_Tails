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

vector <double> centrifugal(double x, double y, double z, double vx, double vy, double vz){
  //function to calculate centrifugal force
  vector <double> omxr, finalx, centri_vector;
  double ang_vel;
	centri_vector.clear();

  ang_vel = omega(m_planet, 1.0);
	omxr = cross_product(0.0, 0.0, ang_vel, x, y, z);

  finalx = cross_product(0.0, 0.0, ang_vel, omxr[0], omxr[1], omxr[2]);

	centri_vector.push_back(finalx[0]);
	centri_vector.push_back(finalx[1]);
	centri_vector.push_back(finalx[2]);

	return centri_vector;
}

vector <double> coriolis(double x, double y, double z, double vx, double vy, double vz){
  //function to calculate coriolis force
  vector <double> omxv,  rvxr, finalx, coriol_vector;

  double ang_vel;

	coriol_vector.clear();

  ang_vel = omega(m_planet, 1.0);
  omxv = cross_product(0.0, 0.0, ang_vel, vx, vy, vz);

	coriol_vector.push_back(omxv[0]);
	coriol_vector.push_back(omxv[1]);
	coriol_vector.push_back(omxv[2]);

	return coriol_vector;

}

vector <double> rad_pressure(double x, double y, double z, double vx, double vy, double vz){
  //function to calculate radiation pressure force, accounting for red shift
  double beta, k, Lum, rad_x, rad_y, rad_z, constant;
  vector <double> rad_vector;

	rad_vector.clear();

	k = opacity(Q_factor);
	Lum = luminosity(Rstar);

	beta = beta_fn(k, Lum, 0.8);


  constant = (beta*G_dim)/(pow(scalar(x- star_x, y, z), 3.0));


	rad_x = constant*(1-((radial_vel(to_vector(vx,vy,vz), sunit_vector(x,y,z)))/c_dim))\
	 *(x-star_x);

	rad_y = constant*(1-((radial_vel(to_vector(vx,vy,vz), sunit_vector(x,y,z)))/c_dim))\
	 *(y);

	rad_z = constant*(1-((radial_vel(to_vector(vx,vy,vz), sunit_vector(x,y,z)))/c_dim))\
	 *(z);


	rad_vector.push_back(rad_x);
	rad_vector.push_back(rad_y);
	rad_vector.push_back(rad_z);

  return rad_vector;

}

vector <double> pr_drag(double x, double y, double z, double vx, double vy, double vz){
  //function to evaluate poynting roberston drag

	vector <double> v_drag, pr_vector;
    pr_vector.clear();
	double constant, beta, k, Lum;

    k = opacity(Q_factor);
	Lum = luminosity(Rstar);


    beta = beta_fn(k, Lum, 0.8);

    constant = (beta*G_dim)/(pow(scalar(x-star_x, y, z), 3.0)*c_dim);

	v_drag = drag_vel(x,y,z,vx,vy,vz);

	pr_vector.push_back(constant* v_drag[0]);
	pr_vector.push_back(constant* v_drag[1]);
	pr_vector.push_back(constant* v_drag[2]);

	return pr_vector;


}
