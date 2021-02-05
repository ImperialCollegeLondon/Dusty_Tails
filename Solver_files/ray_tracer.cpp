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


      for( unsigned int i = 1; i <= cell_no; i++){
          r_a [i] = r_a [i-1] + dr;
          theta_a [i] = theta_a [i-1] + dtheta;
          phi_a [i] = phi_a [i-1] + dphi;

      }

      for (unsigned int i = 0; i <= (cell_no -1); i++){
          r_b [i] = r_a [i] + dr/2.;
          theta_b [i] = theta_a[i] + dtheta/2.;
          phi_b [i] = phi_a [i] + dphi/2.;
      }



}



void calculation_ext(vector <Particle>& particles, double ext [cell_no][cell_no][cell_no], \
                    int nparticles [cell_no][cell_no][cell_no]){
        vector <double> sphere_pos(3, 0.0);
        vector <double> scaled_pos(3, 0.0);
        int r_it, theta_it, phi_it;

        double old_ext;
        double op;
        double pib = 0.;

        double n_mini = 1.0e+23;

        vector <int> r_index, theta_index, phi_index;

        vector <double> p_rs, p_thetas, p_phis;

        vector <double> r_deltas, theta_deltas, phi_deltas;

        double vol_element, partial_vol;

        cout << "total no particles " << particles.size() << endl;

        for( Particle& p : particles) {

              //cout << "in calculation test " << endl;

              sphere_pos = to_spherical(p.position[0], p.position[1], p.position[2]);

              scaled_pos = grid_scaling(sphere_pos);

              if (scaled_pos[0] > 0.0){

              r_it = floor(scaled_pos[0]);
              theta_it = floor(scaled_pos[1]);
              phi_it = floor(scaled_pos[2]);

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

                      op = (3./4.)*(1./4000.)*(1./(p.p_size * 1.e-2));

                      nparticles [r_index[i]][theta_index[j]][phi_index[k]] = nparticles [r_index[i]][theta_index[j]][phi_index[k]] + 1;
                      //cout << nparticles [r_index[i]][theta_index[j]][phi_index[k]] << endl;
                      old_ext = ext [r_index[i]][theta_index[j]][phi_index[k]];
                      ext [r_index[i]][theta_index[j]][phi_index[k]] = old_ext + (partial_vol/ vol_element) * ((n_mini * p.p_mass * 1.0e-3 * op) / (vol_element * pow(a, 3.)));


                  }
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


//verify if particle is within space where optical depth is relevant
//CHANGE THE LIMITS TO THE LIMITS IN EACH VARIABLE
  if ((r<0.) || (r>(n_cells - 1.))){
    scaled = {-1., -1., -1.};

  } else if ((theta<0.) || (theta>(n_cells - 1.))){
    scaled = {-1., -1., -1.};

  } else if ((phi<0.) || (phi>(n_cells - 1.))){
    scaled = {-1., -1., -1.};
  } else {
    scaled = {r, theta, phi};
  }

  return scaled;
}

void optical_depth_test(double ext [cell_no][cell_no][cell_no], double od [cell_no][cell_no][cell_no]){
    for (unsigned int k = 1; k < cell_no; k++){
        for (unsigned int j = 1; j < cell_no; j++){
            for (unsigned int  i = 1; i < cell_no; i++){
                  od[i][j][k] = od[i-1][j][k] + ext[i-1][j][k] * r_reverse(d_dr) * a ;
                }
                }
            }

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
