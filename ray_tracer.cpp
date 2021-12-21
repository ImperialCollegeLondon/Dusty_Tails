#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>
#include <tgmath.h>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include <random>
#include <chrono>
#include "spline.h"

using namespace std;


void build_grids(double *r_a, double *r_b, double *theta_a, \
                double *theta_b, double dr, double dtheta, double dphi,
                double *phi_a, double *phi_b, double r_start, double theta_start, double phi_start){

      r_a[0] = r_start;
      theta_a[0] = theta_start;
      phi_a[0] = phi_start;


      for( unsigned int i = 1; i <= r_cells; i++){
          r_a [i] = r_a [i-1] + dr;
      }
      for( unsigned int i = 1; i <= t_cells; i++){
          theta_a [i] = theta_a [i-1] + dtheta;
      }
      for( unsigned int i = 1; i <= p_cells; i++){
          phi_a [i] = phi_a [i-1] + dphi;
      }


      for (unsigned int i = 0; i <= (r_cells -1); i++){
          r_b [i] = r_a [i] + dr/2.;
      }
      for (unsigned int i = 0; i <= (t_cells -1); i++){
          theta_b [i] = theta_a[i] + dtheta/2.;
      }
      for (unsigned int i = 0; i <= (p_cells -1); i++){
          phi_b [i] = phi_a [i] + dphi/2.;
      }



}


void calculation_ext(vector <Particle>& particles, double ext [r_cells][t_cells][p_cells], double delta_t){
        vector <double> sphere_pos(3, 0.0);
        vector <double> scaled_pos(3, 0.0);
        int r_it, theta_it, phi_it;

        double old_ext;
        double op;
        double pib = 0.;

        //double n_mini = 4.0e+23;
        double mbig = (mdot * T * 0.01) / 1000. ; // 0.01 dependent on when particles are being thrown out of planet
        double n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(1.0e-4, 3));
        //cout <<  "Number of mini particles: " <<  n_mini << endl;

        vector <int> r_index, theta_index, phi_index;

        vector <double> p_rs, p_thetas, p_phis;

        vector <double> r_deltas, theta_deltas, phi_deltas;

        double vol_element, partial_vol;

        for( Particle& p : particles) {

              //cout << "in calculation test " << endl;

              sphere_pos = pos_to_spherical(p.position[0], p.position[1], p.position[2]);

              scaled_pos = grid_scaling(sphere_pos);

              if (scaled_pos[0] > 0.0){

              r_it = floor(scaled_pos[0]);
              //cout << "r_it " << r_it << endl;
              theta_it = floor(scaled_pos[1]);
              //cout << "theta_it " << theta_it << endl;
              phi_it = floor(scaled_pos[2]);
              //cout << "phi_it " << phi_it << endl;

              //particles volume positions in "real" values

              p_rs = { r_reverse(scaled_pos[0] - (dr/2.)), r_reverse(scaled_pos[0] + (dr/2.))};
              p_thetas = { theta_reverse(scaled_pos[1] - (dtheta/2.)), theta_reverse(scaled_pos[1] + (dtheta/2.))};
              p_phis= { phi_reverse(scaled_pos[2] - (dphi/2.)), phi_reverse(scaled_pos[2] + (dphi/2.))};

              if (scaled_pos[0] > r_b[r_it]) {
                r_index = {r_it, r_it + 1};
                r_deltas = {pow(r_reverse(r_a [r_it +1]), 3.0) - pow(p_rs[0], 3.), pow(p_rs[1], 3.) - pow(r_reverse(r_a [r_it +1]), 3.) };
              } else {
                r_index = {r_it - 1, r_it};
                r_deltas = {pow(r_reverse(r_a [r_it]), 3.) - pow(p_rs[0], 3.0) , pow(p_rs[1], 3.0) - pow(r_reverse(r_a [r_it]), 3.) };

              }

              if (scaled_pos[1] > theta_b[theta_it]){
                theta_index = {theta_it, theta_it + 1};
                theta_deltas = {(-cos(theta_reverse(theta_a [theta_it + 1])) + cos(p_thetas[0])), -cos(p_thetas[1]) + cos(theta_reverse(theta_a [theta_it +1])) };
              } else {
                theta_index = {theta_it - 1, theta_it};
                theta_deltas = {(-cos(theta_reverse(theta_a [theta_it])) + cos(p_thetas[0])), -cos(p_thetas[1]) + cos(theta_reverse(theta_a [theta_it]))};
              }

              if (scaled_pos[2] > phi_b[phi_it]){
                phi_index = {phi_it, phi_it + 1};
                phi_deltas = {(phi_reverse(phi_a [phi_it +1]) - p_phis[0]), p_phis[1] - phi_reverse(phi_a [phi_it +1]) };
              } else {
                phi_index = {phi_it - 1, phi_it};
                phi_deltas = {(phi_reverse(phi_a [phi_it]) - p_phis[0]), p_phis[1] - phi_reverse(phi_a [phi_it]) };
              }

              for (unsigned int i = 0; i < 2; i++){
                for (unsigned int j = 0; j < 2; j++){
                  for (unsigned int k = 0; k < 2; k++){

                      //cout << "partial vol: " << endl;
                      partial_vol = 1./3. * abs(r_deltas[i] * theta_deltas[j] * phi_deltas[k]);
                      //cout << partial_vol << endl;
                      vol_element = 1./3. * abs((pow(p_rs[1], 3.) - pow(p_rs[0], 3.)) * (-cos(p_thetas[1]) + cos(p_thetas[0])) * (p_phis[1] - p_phis[0]));
                      //cout << "vol element: " << endl;
                      //cout << vol_element << endl;
                      //op = (3./4.)*(1./4000.)*(1./(p.p_size * 1.e-2)); //opacity in SI units m2/kg
                      op = opacity(p.p_size, p.position[0], p.position[1], p.position[2]); //opacity in CGS units cm2/g
                      op = op*(pow(10.,-4)/pow(10.,-3)); //opacity in SI units
                      //cout << nparticles [r_index[i]][theta_index[j]][phi_index[k]] << endl;
                      old_ext = ext [r_index[i]][theta_index[j]][phi_index[k]];
                      ext [r_index[i]][theta_index[j]][phi_index[k]] = old_ext + (partial_vol/ vol_element) * ((n_mini * p.p_mass*pow(10.,-3) * op) / (vol_element * pow(a, 3)));


                  }
                }
              }
        }
      }


}

vector <double> grid_scaling(vector <double> s_position){
  double r, theta, phi, delta_r, delta_t, delta_p;
  vector <double> scaled(3, 0.0);

  delta_r = d_r_max - d_r_min;
  delta_t = d_t_max - d_t_min;
  delta_p = d_p_max - d_p_min;

  r = (r_cells_d/ (d_r_max - d_r_min))* s_position[0]  - (r_cells_d / (d_r_max - d_r_min)) * d_r_min;
  theta = (t_cells_d/ (d_t_max - d_t_min))* s_position[1]  - (t_cells_d / (d_t_max - d_t_min)) * d_t_min;
  phi = (p_cells_d/ (d_p_max - d_p_min))* s_position[2]  - (p_cells_d /(d_p_max - d_p_min)) * d_p_min;


//verify if particle is within space where optical depth is relevant
//CHANGE THE LIMITS TO THE LIMITS IN EACH VARIABLE
  if ((r<0.) || (r>(r_cells_d - 1.))){
    scaled = {-1., -1., -1.};

  } else if ((theta<0.) || (theta>(t_cells_d - 1.))){
    scaled = {-1., -1., -1.};

  } else if ((phi<0.) || (phi>(p_cells_d - 1.))){
    scaled = {-1., -1., -1.};
  } else {
    scaled = {r, theta, phi};
  }
  return scaled;
}



void optical_depth_calc(double ext [r_cells][t_cells][p_cells], double od [r_cells][t_cells][p_cells]){
    for (unsigned int i = 1; i < r_cells; i++){
        for (unsigned int j = 1; j < t_cells; j++){
            for (unsigned int  k = 1; k < p_cells; k++){
              if (od[i-1][j][k] + ext[i-1][j][k] * d_dr * a > 1.0e-300) {
                  od[i][j][k] = od[i-1][j][k] + ext[i-1][j][k] * d_dr * a ;}
                }
                }
            }
}
//reverse function are from grid scale to "real" scale (in terms of semimajor)
double r_reverse(double old_r){
  double new_r, delta_r;
  delta_r = d_r_max - d_r_min;
  new_r = (old_r/ (r_cells_d/delta_r) ) + d_r_min;
  return new_r;
}

double theta_reverse(double old_theta){
  double new_theta, delta_t;
  delta_t = d_t_max - d_t_min;
  new_theta = (old_theta / (t_cells_d/delta_t)) + d_t_min;
  return new_theta;
}

double phi_reverse(double old_phi){
  double new_phi, delta_p;
  delta_p = d_p_max - d_p_min;
  new_phi = (old_phi / (p_cells_d/delta_p)) + d_p_min;
  return new_phi;
}

vector < vector < vector <double> > >  tau_to_vector(double tau[r_cells][t_cells][p_cells]) {
  //cout << "inside tau to vector" << endl;
  vector < vector < vector <double> > > tauv ;
  for (unsigned int i = 0; i < r_cells; i++){
        tauv.push_back({});
        for (unsigned int j = 0; j < t_cells; j++){
            tauv[i].push_back({});
            for (unsigned int  k = 0; k < p_cells; k++){
              tauv[i][j].push_back(tau[i][j][k]);
              //cout << tau[i][j][k] << endl;
              //cout << tauv[i][j][k] << endl;
            }
        }
  } 
  return tauv;
}

vector <double> r_grid_to_vector(double r[r_cells+1]){
  vector <double> r_v;
  for (unsigned int i = 0; i <r_cells; i++){
    r_v.push_back(r[i]);
  }
  return r_v;
}

vector <double> t_grid_to_vector(double t[t_cells+1]){
  vector <double> t_v;
  for (unsigned int i = 0; i <t_cells; i++){
    t_v.push_back(t[i]);
  }
  return t_v;
}

vector <double> p_grid_to_vector(double p[p_cells+1]){
  vector <double> p_v;
  for (unsigned int i = 0; i <p_cells; i++){
    p_v.push_back(p[i]);
  }
  return p_v;
}

vector <double> vel_grid_scaling(vector <double> s_velocity){
  double vr, vtheta, vphi, delta_r, delta_t, delta_p;
  vector <double> scaled(3, 0.0);

  delta_r = d_r_max - d_r_min;
  delta_t = d_t_max - d_t_min;
  delta_p = d_p_max - d_p_min;

  vr = (r_cells_d/ (d_r_max - d_r_min))* s_velocity[0];
  vtheta = (t_cells_d/ (d_t_max - d_t_min))* s_velocity[1];
  vphi = (p_cells_d/ (d_p_max - d_p_min))* s_velocity[2];

  scaled = {vr, vtheta, vphi};
  return scaled;
}