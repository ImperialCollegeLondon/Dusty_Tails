#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>

using namespace std;

ofstream ofile("ray_tracer_test.bin", ios::out | ios::binary);

#define PI 3.14159

void build_grids(double *r_a_grid, double *r_b_grid, double *theta_a_grid, \
                double *theta_b_grid, double dr, double dtheta, double dphi,
                double *phi_a_grid, double *phi_b_grid, double r_start, double theta_start, double phi_start);

double gauss(double var, double var_c, double std);

void extinction_test(double r_a, double r_b, double theta_a, double theta_b, double phi_a, double phi_b, double function[100][100][100]);

void optical_depth_test(double *extinction, double od[101][101][101]);


double r_a [101];
double r_b [100];
double theta_a [101];
double theta_b [100];
double phi_a [101];
double phi_b [100];

double extinction [100][100][100];
double optical_depth [101][101][101];


double n_cells = 100.0;
double r_min = 0.0;
double r_max = 1.0;
double theta_min = 0;
double theta_max =  1.0*PI;
double phi_min = 0.0;
double phi_max = 2.0*PI;

double dr = ((r_max - r_min)/ n_cells);
double dtheta = (theta_max - theta_min ) / n_cells;
double dphi = ( phi_max - phi_min) / n_cells;

void build_grids(double *r_a, double *r_b, double *theta_a, \
                double *theta_b, double dr, double dtheta, double dphi,
                double *phi_a, double *phi_b, double r_start, double theta_start, double phi_start){

      r_a[0] = r_start;
      theta_a[0] = theta_start;
      phi_a[0] = phi_start;


      for( unsigned int i = 1; i <= 100; i++){
          r_a [i] = r_a [i-1] + dr;
          theta_a [i] = theta_a [i-1] + dtheta;
          phi_a [i] = phi_a [i-1] + dphi;

      }

      for (unsigned int i = 0; i <= 99; i++){
          r_b [i] = r_a [i] + dr/2.;
          theta_b [i] = theta_a[i] + dtheta/2.;
          phi_b [i] = phi_a [i] + dphi/2.;
      }



}

void extinction_test(double *r_a, double *r_b, double *theta_a, double *theta_b, double *phi_a, double *phi_b, double function[100][100][100]){
          double r_gauss, theta_gauss, phi_gauss;
          for (unsigned int i = 0; i< 100; i++){
                for (unsigned int j = 0; j <100; j++){
                    for (unsigned int k = 0; k < 100; k++){
                        r_gauss = gauss(r_b [i], 0.5, 0.05);
                        theta_gauss = gauss(theta_b [j], PI/2.0, 0.5);
                        phi_gauss = gauss(phi_b [k], PI, 0.5);

                        function [i][j][k] = r_gauss * theta_gauss * phi_gauss;

                    }
                }
          }
}

void optical_depth_test(double ext[100][100][100], double od[101][101][101]){
  for (unsigned int i = 1; i <= 100; i++){
        for (unsigned int j = 1; j <= 100; j++){
            for (unsigned int  k = 1; k <= 100; k++){

                  od[i][j][k] = od[i-1][j][k] + ext[i-1][j][k] *dr;
                }


            }
        }
  }



double gauss(double var, double var_c, double std){
      return exp(-pow(var - var_c, 2.0) / (2.0 * pow(std, 2.0)));
}


int main(){
              build_grids(r_a, r_b, theta_a, theta_b, dr, dtheta, dphi, phi_a, phi_b, r_min, theta_min, phi_min);
              extinction_test(r_a, r_b, theta_a, theta_b, phi_a, phi_b, extinction);
              optical_depth_test(extinction, optical_depth);


              for (unsigned int i = 0; i< 101; i++){
                    for (unsigned int j = 0; j <101; j++){
                        for (unsigned int k = 0; k < 101; k++){

                            //cout << extinction [i][j][k] << endl;
                            ofile.write((char*) &r_a[i], sizeof(double));
                            ofile.write((char*) &theta_a[i], sizeof(double));
                            ofile.write((char*) &phi_a[i], sizeof(double));
                            ofile.write((char*) &optical_depth [i][j][k], sizeof(double));


                        }
                    }
              }



}
