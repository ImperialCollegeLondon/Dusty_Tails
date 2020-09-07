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

using namespace std;


void build_grids(double *r_a, double *r_b, double *theta_a, \
                double *theta_b, double dr, double dtheta, double dphi,
                double *phi_a, double *phi_b, double r_start, double theta_start, double phi_start){

      r_a[0] = r_start;
      theta_a[0] = theta_start;
      phi_a[0] = phi_start;


      for( unsigned int i = 1; i <= 200; i++){
          r_a [i] = r_a [i-1] + dr;
          theta_a [i] = theta_a [i-1] + dtheta;
          phi_a [i] = phi_a [i-1] + dphi;

      }

      for (unsigned int i = 0; i <= 199; i++){
          r_b [i] = r_a [i] + dr/2.;
          theta_b [i] = theta_a[i] + dtheta/2.;
          phi_b [i] = phi_a [i] + dphi/2.;
      }



}



void calculation_ext(vector <Particle>& particles, double ext [200][200][200]){
        vector <double> sphere_pos(3, 0.0);
        vector <double> scaled_pos(3, 0.0);
        int r_it, theta_it, phi_it;

        double old_ext;
        double op;

        double n_mini = 1.0e+24;

        vector <int> r_index, theta_index, phi_index;

        vector <double> p_rs, p_thetas, p_phis;

        vector <double> r_deltas, theta_deltas, phi_deltas;

        double vol_element, partial_vol;

        for( Particle& p : particles) {

              cout << "in calculation test " << endl;

              sphere_pos = to_spherical(p.position[0], p.position[1], p.position[2]);
              scaled_pos = grid_scaling(sphere_pos);

              r_it = floor(scaled_pos[0]);
              theta_it = floor(scaled_pos[1]);
              phi_it = floor(scaled_pos[2]);


              //particles volume positions in "real" values

              p_rs = { r_reverse(scaled_pos[0] - (dr/2.)), r_reverse(scaled_pos[0] + (dr/2.))};
              cout << "p_rs " << p_rs[0] << "   " << p_rs[1] << endl;
              p_thetas = { theta_reverse(scaled_pos[1] - (dtheta/2.)), theta_reverse(scaled_pos[1] + (dtheta/2.))};
              cout << "p_thetas " << p_thetas[0] << "   " << p_thetas[1] << endl;
              p_phis= { phi_reverse(scaled_pos[2] - (dphi/2.)), phi_reverse(scaled_pos[2] + (dphi/2.))};
              cout << "p_phis " << p_phis[0] << "   " << p_phis[1] << endl;


              if (scaled_pos[0] > r_b[r_it]) {
                r_index = {r_it, r_it + 1};
                cout << "r_a + 1 " << r_reverse(r_a[r_it +1]) << endl;
                r_deltas = {pow(r_reverse(r_a [r_it +1]), 3.0) - pow(p_rs[0], 3.), pow(p_rs[1], 3.) - pow(r_reverse(r_a [r_it +1]), 3.) };
              } else {
                r_index = {r_it - 1, r_it};
                cout << "r_a " << r_reverse(r_a[r_it]) << endl;
                r_deltas = {pow(r_reverse(r_a [r_it]), 3.) - pow(p_rs[0], 3.0) , pow(p_rs[1], 3.0) - pow(r_reverse(r_a [r_it]), 3.) };

              }

              if (scaled_pos[1] > theta_b[theta_it]){
                theta_index = {theta_it, theta_it + 1};
                theta_deltas = {(-cos(theta_reverse(r_a [theta_it + 1])) + cos(p_thetas[0])), -cos(p_thetas[1]) + cos(theta_reverse(r_a [theta_it +1])) };
              } else {
                theta_index = {theta_it - 1, theta_it};
                theta_deltas = {(-cos(theta_reverse(r_a [theta_it])) + cos(p_thetas[0])), -cos(p_thetas[1]) + cos(theta_reverse(r_a [theta_it]))};
              }

              if (scaled_pos[2] > phi_b[phi_it]){
                phi_index = {phi_it, phi_it + 1};
                phi_deltas = {(phi_reverse(r_a [phi_it +1]) - p_phis[0]), p_phis[1] - phi_reverse(r_a [phi_it +1]) };
              } else {
                phi_index = {phi_it - 1, phi_it};
                phi_deltas = {(phi_reverse(r_a [phi_it]) - p_phis[0]), p_phis[1] - phi_reverse(r_a [phi_it]) };
              }

              for (unsigned int i = 0; i < 2; i++){
                for (unsigned int j = 0; j < 2; j++){
                  for (unsigned int k = 0; k < 2; k++){

                      partial_vol = 1./3. * abs(r_deltas[i] * theta_deltas[j] * phi_deltas[k]);

                      vol_element = 1./3. * abs((pow(p_rs[1], 3.) - pow(p_rs[0], 3.)) * (-cos(p_thetas[1]) + cos(p_thetas[0])) * (p_phis[1] - p_phis[0]));

                      op = (3./4.)*(1./4000)*(1./(p.p_size * 1e-2));

                      cout << "partial vol " << partial_vol << endl;
                      cout << "vol element" << vol_element << endl;

                      cout << "vol fraction " << partial_vol/vol_element << endl;
                      cout << "opacity " << op << endl;

                      old_ext = ext [r_index[i]][theta_index[j]][phi_index[k]];
                      ext [r_index[i]][theta_index[j]][phi_index[k]] = old_ext + (partial_vol/ vol_element) * ((n_mini * p.p_mass * 1.0e-3 * op) / (vol_element * pow(a, 2.)));
                  }
                }
              }
        }

}

vector <double> grid_scaling(vector <double> s_position){
  double r, theta, phi;
  vector <double> scaled(3, 0.0);

  r = (n_cells/ (d_r_max - d_r_min))* s_position[0]  - (n_cells / (d_r_max - d_r_min)) * d_r_min;
  theta = (n_cells/ (d_t_max - d_t_min))* s_position[1]  - (n_cells / (d_t_max - d_t_min)) * d_t_min;
  phi = (n_cells/ (d_p_max - d_p_min))* s_position[2]  - (n_cells /(d_p_max - d_p_min)) * d_p_min;

  scaled = {r, theta, phi};

  return scaled;
}

double r_reverse(double old_r){
  double new_r;

  new_r = (old_r + ((n_cells / (d_r_max - d_r_min)) * d_r_min)) / ((n_cells/ (d_r_max - d_r_min)));

  return new_r;
}

double theta_reverse(double old_theta){
  double new_theta;

  new_theta = (old_theta + ((n_cells / (d_t_max - d_t_min)) * d_t_min)) / ((n_cells/ (d_t_max - d_t_min)));

  return new_theta;
}

double phi_reverse(double old_phi){
  double new_phi;

  new_phi = (old_phi + ((n_cells / (d_p_max - d_p_min)) * d_p_min)) / ((n_cells/ (d_p_max - d_p_min)));

  return new_phi;
}


void extinction_test(double *r_a, double *r_b, double *theta_a, double *theta_b, double *phi_a, double *phi_b, double function[200][200][200]){
          double r_gauss, theta_gauss, phi_gauss;
          for (unsigned int i = 0; i< 200; i++){
                for (unsigned int j = 0; j <200; j++){
                    for (unsigned int k = 0; k < 200; k++){
                        r_gauss = gauss(r_b [i], 0.5, 0.1);
                        theta_gauss = gauss(theta_b [j], PI/2.0, 0.1);
                        phi_gauss = gauss(phi_b [k], PI, 0.1);

                        function [i][j][k] = r_gauss * theta_gauss * phi_gauss;

                    }
                }
          }
}

void od_analytic(double ods[200][200][200]){
  double r_gauss, theta_gauss, phi_gauss;
  for (unsigned int k = 0; k< 200; k++){
        for (unsigned int j = 0; j <200; j++){
            for (unsigned int i = 0; i < 200; i++){
                theta_gauss = gauss(theta_b [j], PI/2.0, 0.1);
                phi_gauss = gauss(phi_b [k], PI, 0.1);

                ods[i][j][k] = theta_gauss * phi_gauss * (-0.1) * pow(PI/2., 1./2.) * erf((0.5-r_b[i])/(pow(2, 1./2.) * 0.1)) + theta_gauss * phi_gauss * (0.1) * pow(PI/2., 1./2.) ;

              }
        }
    }
}

void optical_depth_test(double ext [200][200][200], double od [200][200][200]){
    for (unsigned int k = 1; k < 200; k++){
        for (unsigned int j = 1; j < 200; j++){
            for (unsigned int  i = 1; i < 200; i++){
                  od[i][j][k] = od[i-1][j][k] + ext[i-1][j][k] * r_reverse(r_a[i]-r_a[i-1]);
                }
                }
            }

  }



void error(double analytic [200][200][200], double numerical [200][200][200], double errors [200][200][200]){
  for (unsigned int k = 0; k < 200; k++){
      for (unsigned int j = 0; j < 200; j++){
          for (unsigned int  i = 0; i < 200; i++){

                    errors [i][j][k] = exp(-1.0* (fabs(analytic[i][j][k] - numerical[i][j][k]) / analytic[i][j][k]) );

              }
          }
      }

}


double gauss(double var, double var_c, double std){
      return exp(-1.0 * pow(var - var_c, 2.0) / (2.0 * pow(std, 2.0)));
}
