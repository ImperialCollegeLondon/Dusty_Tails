#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include "spline.h"
#include <random>
#include <chrono>
#include <stdlib.h>
#include <cstring>
#include "opacities.h"
#include <omp.h>
#include <iomanip>

using namespace std;

ofstream output_lc("./data/light_curve.bin", ios::out | ios::binary);

int h_cells, v_cells;
double z_min, z_max, y_min, y_max;
double dh, dv; 

vector< double > h_grid;
vector<double> v_grid;
vector<vector<double>> patches;
vector<vector<double>> taus;

void build_lc_grid(vector<double> &h_grid, vector<double> &v_grid, int h_cells, int v_cells, 
                  double y_min, double z_min, double dh, double dv){


   h_grid[0] = y_min;
   for (unsigned int i=1;i<=h_cells; i++){
      h_grid[i] = h_grid[i-1] + dh;  

     
   }
   v_grid[0] = z_min;
   for (unsigned int k=1;k<=v_cells; k++){
      v_grid[k] = v_grid[k-1] + dv;
   }

}

vector <double> scaled_pos_lc(double x, double y, int h_cells, int v_cells, 
double z_max, double z_min, double y_max, double y_min){
   double x_scaled, y_scaled;
   x_scaled = (h_cells/(y_max-y_min))*x - ((h_cells*y_min)/(y_max-y_min));
   y_scaled = (v_cells/(z_max-z_min))*y - ((v_cells*z_min)/(z_max-z_min));

   return {x_scaled, y_scaled};
}

void  extinction_lc( vector <Particle>& particles, vector <vector <double>> &patches, 
         vector<double> &h_grid, vector<double> &v_grid,
         vector<vector <double>> &taus, int h_cells, int v_cells,
         double z_max, double z_min, double y_max, double y_min, double r_star_a){
         
         
         for ( Particle& p : particles) {
         double rp, hp, vp, dA, dA_cell, cross_sec, shadow;
         vector <double> spos;
         vector <double> h_deltas(2);
         vector <double> v_deltas(2);
         vector <int> h_index, v_index;
         
         int hit, vit;
         if (p.pos_dp[0] > 0.0) {
         
         hp = p.pos_dp[1];
         vp = p.pos_dp[2];
         rp = pow(pow(hp, 2.0) + pow(vp, 2.0), 0.5);
         
         if (rp < r_star_a) {
            spos = scaled_pos_lc(hp, vp, h_cells, v_cells, z_max, z_min, y_max, y_min);
            
            hit = floor(spos[0]);
            vit = floor(spos[1]);
            
            if (hit < h_cells && hit >0  ) {
               if (vit < v_cells && vit >0  ) {
               
            
            
            if ((hit == 0) && (hp-(dh/2.) < h_grid[0])){
                  h_index = {-1, 0};
                  h_deltas[0] = 0.0;
                  h_deltas[1] = abs(abs(hp+(dh/2.)) -abs(h_grid[0]));
            }

             else if ((hit == (h_cells -1)) && (hp+(dh/2.) > h_grid[h_cells])) {
               h_index = {h_cells, -1};
               h_deltas[0] = abs( abs(h_grid[h_cells]) - abs(hp-(dh/2.)));
               h_deltas[1] = 0.0;
            }

             else if (hp > (h_grid[hit]+(dh/2.)) ) {
                h_index = {hit, hit+1};
                h_deltas[0] = abs( abs(h_grid[hit+1]) - abs(hp-(dh/2.)));
                h_deltas[1] = abs( abs(hp+(dh/2.)) - abs(h_grid[hit+1]));

             } 

             else if (hp < (h_grid[hit]+(dh/2.))) {
                
               h_index = {hit-1, hit};
               h_deltas[0] = abs( abs(h_grid[hit]) - abs(hp - (dh/2.)));
               h_deltas[1] = abs( abs(hp+(dh/2.)) - abs(h_grid[hit]));

             }
            
            if ((vit == 0) && (vp-(dv/2.) < v_grid[0])){
                  v_index = {-1, 0};
                  v_deltas[0] = 0.0;
                  v_deltas[1] = abs(abs(vp+(dv/2.)) -abs(v_grid[0]));
            }
            else if ((vit == (v_cells -1)) && (vp+(dv/2.) > v_grid[v_cells])) {
               v_index = {v_cells, -1};
               v_deltas[0] = abs( abs(v_grid[v_cells]) - abs(vp-(dv/2.)));
               v_deltas[1] = 0.0;
            }

            else if (vp > (v_grid[vit]+(dv/2.)) ) {
               v_index = {vit, vit+1};
               v_deltas[0] = abs( abs(v_grid[vit+1]) - abs(vp-(dv/2.)));
               v_deltas[1] = abs( abs(vp+(dv/2.)) - abs(v_grid[vit+1]));

            } 

            else if (vp < (v_grid[vit]+(dv/2.))) {
               v_index = {vit-1, vit};
               v_deltas[0] = abs( abs(v_grid[vit]) - abs(vp - (dv/2.)));
               v_deltas[1] = abs( abs(vp+(dv/2.)) - abs(v_grid[vit]));

            }
            
            for (unsigned int i=0; i<2; i++){
               for (unsigned int j=0; j<2; j++){
                     
                     dA = h_deltas[i]*v_deltas[j];
                    
                     if (dA > 0.0) {
                     dA_cell = patches[h_index[i]][v_index[j]];
                     cross_sec = (p.opac_planck * p.mass) / (pow(a, 2.0));
                     
                     shadow = ((dA/(dh*dv)) * cross_sec *p.n_mini) / dA_cell;
                     if (shadow < 1.0e-20) {
                        shadow = 0.0;
                     }
                     
                     taus[h_index[i]][v_index[j]] = taus[h_index[i]][v_index[j]] + shadow;
                     
                     }
               }
            }
            } 
            } 
         }
      }  
         }
   }


void grid_cells_lc(vector<double> &h_grid, vector<double> &v_grid, 
               vector<vector<double>> &patches, 
               int h_cells, int v_cells, double r_star_a){
                  
   double p1[2], p2[2], p3[2], p4[2];
   double r1, r2, r3, r4;
   double check;
   double delta_x, delta_y;
   double h, a, b;

   
   for (int m=0; m<h_cells; m++){
      for (int n=0; n<v_cells; n++){
            //cout <<  m << " " << n << " " << endl;
            check = 0;
            p1[0] = h_grid[m];
            //cout << "x1 " << p1[0] << " y1 " << p1[1] << endl;
            p1[1] = v_grid[n];
            p2[0] = h_grid[m+1];

            p2[1] = v_grid[n];
            //cout << "x2 " << p2[0] << " y2 " << p2[1] << endl;
            p3[0] = h_grid[m+1];
            p3[1] = v_grid[n+1];
             //cout << "x3 " << p3[0] << " y3 " << p3[1] << endl;
            p4[0] = h_grid[m];
            p4[1] = v_grid[n+1];
            //cout << "x4 " << p4[0] << " y4 " << p4[1] << endl;
           
           
            r1 = pow(pow(p1[0],2) + pow(p1[1],2), 0.5);
            r2 = pow(pow(p2[0],2) + pow(p2[1],2), 0.5);
            r3 = pow(pow(p3[0],2) + pow(p3[1],2), 0.5);
            r4 = pow(pow(p4[0],2) + pow(p4[1],2), 0.5);

            //cell outside

            
            if ((r1 >=r_star_a) && (r2 >= r_star_a) && (r3 >=r_star_a) && (r4 >= r_star_a)){
               // cout << "cell out" << endl;
               patches[m][n] = 0.0;
            }

            //cell inside
            else if ((r1 < r_star_a) && (r2 < r_star_a) && (r3 < r_star_a) && (r4 < r_star_a)) {
               patches[m][n] = dh*dv;
               
            }

            else if (r1 <r_star_a) {
               //Triangle shape, 1st quadrant, 1 in, 3 out
               if ((r2 >= r_star_a) && (r3>=r_star_a) && (r4 >=r_star_a) ) {
                  delta_x = pow(pow(r_star_a,2)-pow(p1[1],2), 0.5) - abs(p1[0]);
                  delta_y = pow(pow(r_star_a,2)-pow(p1[0],2), 0.5) - abs(p1[1]);
                  patches[m][n] = (delta_x * delta_y)/2.;
                  
               } 
               //1st quadrant, 3 in, 1 out
               else if ((r2 < r_star_a) && (r4 < r_star_a) && (r3 >=r_star_a)) {
                  delta_x = abs(p3[0]) - pow(pow(r_star_a,2)-pow(p3[1],2), 0.5);
                  delta_y = abs(p3[1]) - pow(pow(r_star_a,2)-pow(p3[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                
               }
            
               //2nd quadrant, 3 in, 1 out
               else if ((r2 < r_star_a) && (r3 < r_star_a) && (r4 >=r_star_a)) {
                  delta_x = abs(p4[0]) - pow(pow(r_star_a,2)-pow(p4[1],2), 0.5);
                  delta_y = abs(p4[1]) - pow(pow(r_star_a,2)-pow(p4[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                  
               }
               //4th quadrant, 3 in, 1 out
               else if ((r3 < r_star_a) && (r4 < r_star_a) && (r2 >= r_star_a)) {
                  delta_x = abs(p2[0]) - pow(pow(r_star_a,2)-pow(p2[1],2), 0.5);
                  delta_y = abs(p2[1]) - pow(pow(r_star_a,2)-pow(p2[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                  
               }
               //trapezium, height parallel to x axis, 1 and 2 in
               else if ((r2 < r_star_a) && (r4 >= r_star_a) && (r3 >= r_star_a)) {
                  h = abs(abs(p2[0])-abs(p1[0]));
                  a = pow(pow(r_star_a,2)-pow(p2[0],2), 0.5) - abs(p2[1]);
                  b = pow(pow(r_star_a,2)-pow(p1[0],2), 0.5) - abs(p1[1]);
                  patches[m][n] = 0.5*h*(a+b);
                  
                  
               }
                //trapezium, height parallel to y axis, 1 and 4 in
               else if ((r4 < r_star_a) && (r2 >= r_star_a) && (r3 >= r_star_a)) {
                  h = abs(abs(p4[1])-abs(p1[1]));
                  a = pow(pow(r_star_a,2)-pow(p4[1],2), 0.5) - abs(p4[0]);
                  b = pow(pow(r_star_a,2)-pow(p1[1],2), 0.5) - abs(p1[0]);
                  patches[m][n] = 0.5*h*(a+b);
                  
               }   
            } 

            else if (r3 < r_star_a) {
               //triangle, 3 in
               if ((r1 >=r_star_a) && (r2>= r_star_a) && (r4 >=r_star_a)) {
                  delta_x = pow(pow(r_star_a,2)-pow(p3[1],2), 0.5) - abs(p3[0]);
                  delta_y = pow(pow(r_star_a,2)-pow(p3[0],2), 0.5) - abs(p3[1]);
                  patches[m][n] = delta_x * delta_y * 0.5;
               }
               //triangle, 3 2 and 4 in, 1 out
               else if ((r2 < r_star_a) && (r4 < r_star_a) && (r1 >=r_star_a)) {
                  delta_x = abs(p1[0]) - pow(pow(r_star_a,2)-pow(p1[1],2), 0.5);
                  delta_y = abs(p1[1]) - pow(pow(r_star_a,2)-pow(p1[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
               }

               else if ((r2 < r_star_a) && (r1 >= r_star_a) && (r4 >= r_star_a)) {
                  h = abs(abs(p3[1])-abs(p2[1]));
                  a = pow(pow(r_star_a,2)-pow(p3[1],2), 0.5) - abs(p3[0]);
                  b = pow(pow(r_star_a,2)-pow(p2[1],2), 0.5) - abs(p2[0]);
                  patches[m][n] = 0.5*h*(a+b);
               }
               else if ((r4 <r_star_a) && (r1 >=r_star_a) && (r2> r_star_a)) {
                  h = abs(abs(p3[0])-abs(p4[0]));
                  a = pow(pow(r_star_a,2)-pow(p3[0],2), 0.5) - abs(p3[1]);
                  b = pow(pow(r_star_a,2)-pow(p4[0],2), 0.5) - abs(p4[1]);
                  patches[m][n] = 0.5*h*(a+b);
               }
            }

            else if (r2 <r_star_a) {
               if ((r1 >=r_star_a) && (r3>=r_star_a) && (r3>=r_star_a)){
                  delta_x = pow(pow(r_star_a,2)-pow(p2[1],2), 0.5) - abs(p2[0]);
                  delta_y = pow(pow(r_star_a,2)-pow(p2[0],2), 0.5) - abs(p2[1]);
                  patches[m][n] = delta_x * delta_y * 0.5;
               }
            }

            else if (r4 < r_star_a) {
               if ((r1 >=r_star_a) && (r2>=r_star_a) && (r3>=r_star_a)){
                  delta_x = pow(pow(r_star_a,2)-pow(p4[1],2), 0.5) - abs(p4[0]);
                  delta_y = pow(pow(r_star_a,2)-pow(p4[0],2), 0.5) - abs(p4[1]);
                  patches[m][n] = delta_x * delta_y * 0.5;
               }
            }
            else {
               cout << "something went wrong " << endl;
               abort;
            }
      }
   }
   
}

double forward_scat(double g, double scat_opac, double mass, double x, double y, double z, double s, double tau, double n_mini){
      double phase;
      double d_obs_star = earth_star*pc;
      double phase_angle;
      double r_dust = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      phase_angle = acos((1/r_dust) * (x*d_obs_star - pow(r_dust,2)) 
                     / sqrt(pow(d_obs_star, 2.) - 2.0*d_obs_star*x + pow(r_dust,2.)));
      phase = ((1-pow(g,2))/pow(1+pow(g,2)-2*g*cos(phase_angle), 3./2.));
     
      return  exp(-tau)*phase * scat_opac*mass*n_mini/ (4.0*PI*pow(x,2.)+pow(y,2.)+pow(z,2.));
}

double flux(vector <vector <double>> &taus, vector< vector <double>> &patches, int h_cells, int v_cells, double r_star_a){
   double f, total_grid_area;
   total_grid_area = 0.;
   for (int m=0; m<h_cells; m++){
      for (int n=0; n<v_cells; n++){
         total_grid_area = total_grid_area + patches[m][n];
      }
   }

   f = 0.0;
   for (int m=0; m<h_cells; m++){
      for (int n=0; n<v_cells; n++){
         if (taus[m][n] != 0.0){
            f = f + patches[m][n] *(1.0- exp(-1.0*taus[m][n]));
           
         }
      
         
      }
   }

   return 1.0 - (f/(PI*pow(r_star_a,2)));
}

void light_curve(vector<Particle>& particles, double current_t){
std::cout << std::fixed;
std::cout << std::setprecision(8) << endl;
std::cout << std::scientific << endl;
using namespace std;
double inclination, Rstar_a;

inclination = acos((Rstar*Rsun_cgs*b_p)/a); //in radians

Rstar_a = (Rstar*Rsun_cgs)/a; //stellar radius in terms of the semimajor axis

double f_test_o = 1.0;
double phi;
double forward_flux = 0.0;
phi = 2.0*PI*(current_t - t0);

for( Particle& p : particles) {
      //first projection
        p.pos_p = {p.position[0] * cos(phi) - p.position[1]*sin(phi),
                  p.position[0] * sin(phi) + p.position[1]*cos(phi),
                   p.position[2]};
        
      //second projection
        p.pos_dp = {p.pos_p[0]*sin(inclination) + p.pos_p[2]*cos(inclination),
                   p.pos_p[1],
                  -1.0*p.pos_p[0]*cos(inclination) + p.pos_p[2]*sin(inclination)};
        
      //forward scattering
        if (p.pos_dp[0]>0.0) {
        p.f_scat = forward_scat(p.gsca, 
                 p.opac_scat, p.mass, p.pos_dp[0]*a, p.pos_dp[1]*a, 
                 p.pos_dp[2]*a, p.size, p.tau_d, p.n_mini);
        forward_flux = forward_flux + p.f_scat;
        } else {
         p.f_scat = 0.0;
        }
}

   v_cells = 80;
   double f_ext = 0.;
   double f_total = 0.;
   double xp_planet, yp_planet, zp_planet;
   double xdp_planet, ydp_planet, zdp_planet;
   xp_planet = cos(phi);
   yp_planet = sin(phi);
   zp_planet = 0.0;
   xdp_planet = xp_planet*sin(inclination);
   ydp_planet = yp_planet;
   zdp_planet = -xp_planet*cos(inclination);

   z_min = zdp_planet - 0.04;
   z_max = zdp_planet + 0.04;
   y_min = -Rstar_a;
   y_max = Rstar_a;

   h_cells = round((y_max-y_min) / ((z_max-z_min)/v_cells));
   dh = (y_max-y_min)/h_cells;
   dv = (z_max-z_min)/v_cells;

   for (int j=0; j<=h_cells; j++) {
      h_grid.push_back(0.0);
   }
   for (int j=0; j<=v_cells; j++) {
      v_grid.push_back(0.0);
   }
   for (int m=0; m<h_cells; m++){
      patches.push_back({});
      for (int n=0; n<v_cells; n++){
         patches[m].push_back(0.0);
      }

   }
   
   build_lc_grid(h_grid, v_grid, h_cells, v_cells, y_min, z_min, dh, dv); 
   grid_cells_lc(h_grid, v_grid, patches, h_cells, v_cells, Rstar_a);

   for (int w=0; w<h_cells; w++){
         taus.push_back({});
         for (int y=0; y<v_cells; y++){
            taus[w].push_back(0.0);
         }
   }
         
          
   extinction_lc(particles, patches,  h_grid, v_grid,taus, h_cells, v_cells, 
         z_max, z_min, y_max, y_min, Rstar_a);
        
   output_lc.write((char*) &current_t, sizeof(double));
         
   f_ext = flux(taus, patches, h_cells, v_cells, Rstar_a);

   cout << "Extinction is " << f_ext << endl;
   cout << "Scattering is " << forward_flux << endl;
   output_lc.write((char*) &f_ext, sizeof(double));
   output_lc.write((char*) &forward_flux, sizeof(double));
   f_total = f_ext + forward_flux;
   cout << "Transit depth is " << f_total << endl;
   output_lc.write((char*) &f_total, sizeof(double));

   taus.clear();
   h_grid.clear();
   v_grid.clear();
   patches.clear();

    }