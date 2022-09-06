#include <iostream>
#include<iostream>
#include<fstream>
#include<vector>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <string>
#include "opacities.h"
#include <iomanip>

using namespace std;

#define amu 1.661e-24
#define kb 1.381e-16

#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define Msun_cgs 1.9885e+33 // solar mass
#define Mearth_cgs 5.972e+27 // earth mass grams
#define Rearth_cgs 6.378e+8 // earth radius cm

struct dust_read {
   double timestamp;
   long int id;
   double x_dust;
   double y_dust;
   double z_dust;
   double vx_dust;
   double vy_dust;
   double vz_dust;
   double s_dust;
   double h_dust;
   double m_dust;
   double temp_dust;
   double tau_dust;
   double kappa_dust_abs;
   double kappa_dust_scat;
   double kappa_planck;
};

struct dust {
   double timestamp;
   double phi;
   double x_p;
   double y_p;
   double z_p;
   double x_dp;
   double y_dp;
   double z_dp;
   double m;
   double kappa;
   double kappa_scat;
   double size;
   double tau;
};
template<typename T>
vector<double> linspace(T start_in, T end_in, int num_in){

  vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}


class Particle {
   public:
    long int id; //unique id of the particle
    vector <double> position; //current position of the particle
    vector <double> velocity; //current velocity of the particle
    vector <double> pos_spherical; // current position in spherical coordinates
    vector <double> v_spherical; //curretn velocity in spherical coordinates
    double size; //size of the particle
    double opac_abs;
    double opac_scat; 
    double opac_planck;
    double h_updated; //current optimal time step for particle
    double mass; //mass of particle
    double tau_d; //optical depth
    double temp_d; //particle temperature

};
#define PI 3.14159
const double gyr = pow(10.,9) * 365. * 24. * 60. * 60.;
double mdot_read; //mass loss rate of planet in Earth masses per Gyr

double s_0, rho_d;
int h_cells, v_cells;
double Temp, T , mbig,n_mini,m_star, a_p, inclination, r_star;
double z_min, z_max, y_min, y_max;
double dh, dv; 

const int timesteps = 100;

double angle[timesteps];
string opacity_dir;
string opac_data;
Opacities opac;

ofstream output("./output.bin", ios::out | ios::binary); 

double* angles(double times[timesteps], double t0, double angle[timesteps]){
   for (int i=0; i<=timesteps; i++){
      angle[i] = 2.0*PI * (times[i] - t0);
   }
   return angle;
}

void build_grid(vector<double> &h_grid, vector<double> &v_grid, int h_cells, int v_cells){
   h_grid[0] = y_min;
   for (unsigned int i=1;i<=h_cells; i++){
      h_grid[i] = h_grid[i-1] + dh;  
     
   }
   v_grid[0] = z_min;
   for (unsigned int k=1;k<=v_cells; k++){
      v_grid[k] = v_grid[k-1] + dv;
   }
   //cout << "grid built " << endl;
}

vector <dust_read> read_data(){
  std::fstream input;
  input.open("./input.bin", std::fstream::in | std::fstream::binary);
  input.seekg(0, ios::end);
  int size=input.tellg();
  input.seekg(0, ios::beg);
  long int total = size/sizeof(dust);
  cout << "total " << total << endl;
  vector <dust_read> dust_grains_out;
  for(int i = 0; i < total; i++){
       dust_grains_out.push_back(dust_read());
       input.read((char *) &dust_grains_out[i], sizeof(dust_read));
       
       if (dust_grains_out[i].s_dust < 0.0) {
         dust_grains_out.pop_back();
         
       }
       
       if (dust_grains_out[i].id == 0) {
          //dust_grains_out.erase(dust_grains_out.end());
          dust_grains_out.pop_back();

          break;
          
       }

    }
   input.close();
   cout << "successfully read the data" << endl;
   return dust_grains_out;
}

vector <double> scaled_pos(double x, double y, int h_cells, int v_cells, 
double z_max, double z_min, double y_max, double y_min){
   double x_scaled, y_scaled;
   x_scaled = (h_cells/(y_max-y_min))*x - ((h_cells*y_min)/(y_max-y_min));
   y_scaled = (v_cells/(z_max-z_min))*y - ((v_cells*z_min)/(z_max-z_min));

   return {x_scaled, y_scaled};
}

void  extinction( vector <dust> particles, vector <vector <double>> &patches, 
         vector<double> &h_grid, vector<double> &v_grid,
         vector<vector <double>> &taus, int h_cells, int v_cells,
         double z_max, double z_min, double y_max, double y_min){
         
         
         for ( dust& p : particles) {
         double rp, hp, vp, dA, dA_cell, sigma, shadow;
         vector <double> spos;
         vector <double> h_deltas(2);
         vector <double> v_deltas(2);
         vector <int> h_index, v_index;
         
         int hit, vit;
         if (p.x_dp > 0.0) {
         hp = p.y_dp;
         vp = p.z_dp;
         rp = pow(pow(hp, 2.0) + pow(vp, 2.0), 0.5);
         if (rp < r_star) {
            spos = scaled_pos(hp, vp, h_cells, v_cells, z_max, z_min, y_max, y_min);
            
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
                     sigma = (p.kappa * p.m) / (pow(a_p, 2.0));
                     
                     shadow = ((dA/(dh*dv)) * sigma*n_mini) / dA_cell;
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

         //cout << "obtain extinction grid successfully" << endl;

   }




void grid_cells(vector<double> &h_grid, vector<double> &v_grid, 
               vector<vector<double>> &patches, 
               int h_cells, int v_cells){
                  
   double p1[2], p2[2], p3[2], p4[2];
   double r1, r2, r3, r4;
   double check;
   double delta_x, delta_y;
   double h, a, b;
   //cout << "at grid cells " << endl;
   
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

            
            if ((r1 >=r_star) && (r2 >= r_star) && (r3 >=r_star) && (r4 >= r_star)){
               // cout << "cell out" << endl;
               patches[m][n] = 0.0;
            }

            //cell inside
            else if ((r1 < r_star) && (r2 < r_star) && (r3 < r_star) && (r4 < r_star)) {
               patches[m][n] = dh*dv;
               
            }

            else if (r1 <r_star) {
               //Triangle shape, 1st quadrant, 1 in, 3 out
               if ((r2 >= r_star) && (r3>=r_star) && (r4 >=r_star) ) {
                  delta_x = pow(pow(r_star,2)-pow(p1[1],2), 0.5) - abs(p1[0]);
                  delta_y = pow(pow(r_star,2)-pow(p1[0],2), 0.5) - abs(p1[1]);
                  patches[m][n] = (delta_x * delta_y)/2.;
                  
               } 
               //1st quadrant, 3 in, 1 out
               else if ((r2 < r_star) && (r4 < r_star) && (r3 >=r_star)) {
                  delta_x = abs(p3[0]) - pow(pow(r_star,2)-pow(p3[1],2), 0.5);
                  delta_y = abs(p3[1]) - pow(pow(r_star,2)-pow(p3[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                
               }
            
               //2nd quadrant, 3 in, 1 out
               else if ((r2 < r_star) && (r3 < r_star) && (r4 >=r_star)) {
                  delta_x = abs(p4[0]) - pow(pow(r_star,2)-pow(p4[1],2), 0.5);
                  delta_y = abs(p4[1]) - pow(pow(r_star,2)-pow(p4[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                  
               }
               //4th quadrant, 3 in, 1 out
               else if ((r3 < r_star) && (r4 < r_star) && (r2 >= r_star)) {
                  delta_x = abs(p2[0]) - pow(pow(r_star,2)-pow(p2[1],2), 0.5);
                  delta_y = abs(p2[1]) - pow(pow(r_star,2)-pow(p2[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
                  
               }
               //trapezium, height parallel to x axis, 1 and 2 in
               else if ((r2 < r_star) && (r4 >= r_star) && (r3 >= r_star)) {
                  h = abs(abs(p2[0])-abs(p1[0]));
                  a = pow(pow(r_star,2)-pow(p2[0],2), 0.5) - abs(p2[1]);
                  b = pow(pow(r_star,2)-pow(p1[0],2), 0.5) - abs(p1[1]);
                  patches[m][n] = 0.5*h*(a+b);
                  
                  
               }
                //trapezium, height parallel to y axis, 1 and 4 in
               else if ((r4 < r_star) && (r2 >= r_star) && (r3 >= r_star)) {
                  h = abs(abs(p4[1])-abs(p1[1]));
                  a = pow(pow(r_star,2)-pow(p4[1],2), 0.5) - abs(p4[0]);
                  b = pow(pow(r_star,2)-pow(p1[1],2), 0.5) - abs(p1[0]);
                  patches[m][n] = 0.5*h*(a+b);
                  
               }   
            } 

            else if (r3 < r_star) {
               //triangle, 3 in
               if ((r1 >=r_star) && (r2>= r_star) && (r4 >=r_star)) {
                  delta_x = pow(pow(r_star,2)-pow(p3[1],2), 0.5) - abs(p3[0]);
                  delta_y = pow(pow(r_star,2)-pow(p3[0],2), 0.5) - abs(p3[1]);
                  patches[m][n] = delta_x * delta_y * 0.5;
               }
               //triangle, 3 2 and 4 in, 1 out
               else if ((r2 < r_star) && (r4 < r_star) && (r1 >=r_star)) {
                  delta_x = abs(p1[0]) - pow(pow(r_star,2)-pow(p1[1],2), 0.5);
                  delta_y = abs(p1[1]) - pow(pow(r_star,2)-pow(p1[0],2), 0.5);
                  patches[m][n] = (dh*dv) - 0.5*delta_x * delta_y;
               }

               else if ((r2 < r_star) && (r1 >= r_star) && (r4 >= r_star)) {
                  h = abs(abs(p3[1])-abs(p2[1]));
                  a = pow(pow(r_star,2)-pow(p3[1],2), 0.5) - abs(p3[0]);
                  b = pow(pow(r_star,2)-pow(p2[1],2), 0.5) - abs(p2[0]);
                  patches[m][n] = 0.5*h*(a+b);
               }
               else if ((r4 <r_star) && (r1 >=r_star) && (r2> r_star)) {
                  h = abs(abs(p3[0])-abs(p4[0]));
                  a = pow(pow(r_star,2)-pow(p3[0],2), 0.5) - abs(p3[1]);
                  b = pow(pow(r_star,2)-pow(p4[0],2), 0.5) - abs(p4[1]);
                  patches[m][n] = 0.5*h*(a+b);
               }
            }

            else if (r2 <r_star) {
               if ((r1 >=r_star) && (r3>=r_star) && (r3>=r_star)){
                  delta_x = pow(pow(r_star,2)-pow(p2[1],2), 0.5) - abs(p2[0]);
                  delta_y = pow(pow(r_star,2)-pow(p2[0],2), 0.5) - abs(p2[1]);
                  patches[m][n] = delta_x * delta_y * 0.5;
               }
            }

            else if (r4 < r_star) {
               if ((r1 >=r_star) && (r2>=r_star) && (r3>=r_star)){
                  delta_x = pow(pow(r_star,2)-pow(p4[1],2), 0.5) - abs(p4[0]);
                  delta_y = pow(pow(r_star,2)-pow(p4[0],2), 0.5) - abs(p4[1]);
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

double forward_scat(double g, double scat_opac, double mass, double x, double y, double z, double s, double tau){
      double phase;
      double d_obs_star = 1.911e+21;
      double phase_angle;
      double r_dust = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      phase_angle = acos((1/r_dust) * (x*d_obs_star - pow(r_dust,2)) 
                     / sqrt(pow(d_obs_star, 2.) - 2.0*d_obs_star*x + pow(r_dust,2.)));
      phase = ((1-pow(g,2))/pow(1+pow(g,2)-2*g*cos(phase_angle), 3./2.));
     
      return  exp(-tau)*phase * scat_opac*mass*n_mini/ (4.0*PI*pow(x,2.)+pow(y,2.)+pow(z,2.));
}

double flux(vector <vector <double>> &taus, vector< vector <double>> &patches, int h_cells, int v_cells){
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

   return 1.0 - (f/(PI*pow(r_star,2)));
}

int main(){


std::cout << std::fixed;
std::cout << std::setprecision(8) << endl;
std::cout << std::scientific << endl;
using namespace std;
#include <cstring>
   string planet, composition;
   string opac_data, line;
   double no_p;
   
  
   int in_c = 0;
   fstream input_dust;
	input_dust.open("./light_curve/lightcurve.in", ios::in);

    if (input_dust.is_open()){
    while ( getline (input_dust,line) )
    {
      if (in_c == 0) {
        planet = line.substr(0,10);
        cout << planet << endl;
      }
      if (in_c == 1) {
        composition = line.substr(0,15);
        cout << composition << endl;
      }
      if (in_c == 2) {
         s_0 = stod(line.substr(0,5)) * 1.0e-4;
         cout << "s_0" << s_0 << endl;
      }
      if (in_c == 3) {
         no_p =  stod(line.substr(0,5));
         cout << no_p << endl;
      }
      if (in_c ==4) {
         mdot_read = stod(line.substr(0,5));
      }
      in_c = in_c + 1;
      
    }
    input_dust.close();
  } else {
     cout << "oops" << endl;
  }

   double mdot = mdot_read*Mearth_cgs/gyr;
   if (composition.substr(0,5) == "Al2O3") {
      cout << "Dust is composed of Corundum." << endl;
      opac_data = "corundum_K95";
      double A = 7.74e+4; //clausius claperyon relation
      double Bp = 39.3; //clausius claperyon relation
      rho_d = 4.0; //dust density
      double mu = 101.961;
      double alpha = 0.1;
   } else if (composition.substr(0,6)=="MgSiO3"){
      cout << "Dust is composed of Enstatite." << endl;
      opac_data = "enstatite_J98_J94_D95";
      double A = 6.89e+4;
      double Bp = 37.8;
      rho_d = 3.20;
      double mu = 100.389;
      double alpha = 0.1;

   }  else if (composition.substr(0,8)=="MgFeSiO4") {
  cout << "Dust is composed of Olivine [MgFe]2SiO4" << endl;
  opac_data = "MgFeSiO4_D95";
  //A and B are the same values as for Forsterite (Perez-Becker+ 2013)
  double A = 6.53e+4; 
  double Bp = 34.3;
  rho_d = 3.8;
  double mu = 172.23;
  double alpha = 0.1;

} else if (composition.substr(0,12)=="Mg08Fe12SiO4") {
  cout << "Dust is composed of Olivine (Mg08,Fe12)" << endl;
  opac_data = "Mg08Fe12SiO4_D95";
  double A = 6.53e+4; 
  double Bp = 34.3;
  rho_d = 3.80;
  double mu = 178.538;
  double alpha = 0.1;
}  else{
   cout << "Composition unknown, stopping.";
   abort();
   }
   if (planet.substr(0,9) == "KIC1255_b") {
      Temp = 4550.0;
      T = 15.68*60.*60.; //planetary period in seconds
      mbig = (mdot * T * 0.01) / no_p ; // 0.01 dependent on when particles are being thrown out of planet

      n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
      cout << "n_mini " << n_mini << endl;
      m_star = 0.66; //stellar mass in solar masses
      a_p = pow((G_cgs*m_star*Msun_cgs* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in cgs
      inclination = 1.425; //in radians
      r_star = 0.2415; //stellar radius in terms of the semimajor axis
   } else {
      Temp = 3830.0;
      T = 9.15*60.*60.; //planetary period in seconds
      mbig = (mdot * T * 0.01) / no_p ; // 0.01 dependent on when particles are being thrown out of planet
      n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
      m_star = 0.60; //stellar mass in solar masses
      a_p = pow((G_cgs*m_star*Msun_cgs* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in cgs
      inclination = acos((0.57*Rsun_cgs*0.68)/a_p); //in radians
      cout << "i " << inclination << endl;
      r_star = (0.57*Rsun_cgs)/a_p; //stellar radius in terms of the semimajor axis
      //cout << r_star << endl;
   }
   vector <double> period_steps = {};
   vector <dust_read> particles_read;
   vector <dust> particles;
   particles_read = read_data();
   double t0 = 0.5;
   int counter = 0;
   double f_test_o = 1.0;
   vector <double> timestamps;
   vector <double> transit_depths = {};
   vector <double> grid_cell_size = {};
   int T_int { static_cast<int> (Temp)};
   string T_int_s = to_string(T_int);
   opacity_dir = "./opacs_jankovic/calc_dust_opac/"+opac_data+"/opac_";
   
   opac.read_data((opacity_dir+"temp.dat").c_str(), (opacity_dir+"sizes.dat").c_str(),
        (opacity_dir+"planck_abs_tstar"+T_int_s+".dat").c_str(), 
        (opacity_dir+"planck_sca_tstar"+T_int_s+".dat").c_str(),
        (opacity_dir+"planck_gsc_tstar"+T_int_s+".dat").c_str(),
        (opacity_dir+"planck_abs.dat").c_str(), (opacity_dir+"planck_sca.dat").c_str(),
        true);
   cout << "read opacity data successfully " << endl;
   double g_test;
   counter = 0;

    for( dust_read& p : particles_read) {
        particles.push_back(dust());
      
        particles[counter].timestamp = p.timestamp;
        double key = p.timestamp;
        if (find(timestamps.begin(), timestamps.end(), key) == timestamps.end()) {
           //cout << key << endl;
           timestamps.push_back(key);
       }
        particles[counter].phi = 2.0*PI * (p.timestamp - t0);
        particles[counter].x_p = p.x_dust * cos(particles[counter].phi) - p.y_dust * sin(particles[counter].phi);
        particles[counter].y_p = p.x_dust * sin(particles[counter].phi) + p.y_dust * cos(particles[counter].phi);
        particles[counter].z_p = p.z_dust;
        particles[counter].x_dp = particles[counter].x_p * sin(inclination) + particles[counter].z_p * cos(inclination);
        particles[counter].y_dp = particles[counter].y_p;
        particles[counter].z_dp = -particles[counter].x_p * cos(inclination) + particles[counter].z_p * sin(inclination);
        particles[counter].m = p.m_dust;
        particles[counter].kappa = p.kappa_planck;
        particles[counter].kappa_scat = p.kappa_dust_scat;
        particles[counter].size = p.s_dust;
        particles[counter].tau = p.tau_dust;

        //cout << p.s_dust << endl;
        //cout << p.kappa_planck << endl;
        //cout << forward_scat(opac.stellar_gsc(p.s_dust), p.kappa_dust_scat, p.m_dust) << endl;;
       
        counter = counter + 1;
        }


   vector< double > h_grid;
   vector<double> v_grid;
   vector<vector<double>> patches;
   vector<vector<double>> taus;
   vector <dust> particles_calc;

   for (int i=120; i<121; i++){

      //memset(taus, 0, sizeof(taus));
      v_cells = i;
      double f_test = 0.;
      double f_test_o = 1.0;
      double x_planet, y_planet, z_planet;
      for ( int it=0; it<timestamps.size(); it++) {
      //for ( int it=0; it<1; it++) {
         cout << timestamps[it] << endl;
      
      double theta;
      
      theta = 2.0*PI*(timestamps[it]-t0);
      double xp_planet, yp_planet, zp_planet;
      double xdp_planet, ydp_planet, zdp_planet;
      xp_planet = cos(theta);
      yp_planet = sin(theta);
      zp_planet = 0.0;
      xdp_planet = xp_planet*sin(inclination);
      ydp_planet = yp_planet;
      zdp_planet = -xp_planet*cos(inclination);

      z_min = zdp_planet - 0.04;
      z_max = zdp_planet + 0.04;

      //cout << "z_planet " << z_planet << endl;
      z_min = zdp_planet - 0.04;
      z_max = zdp_planet + 0.04;
      y_min = -r_star;

      y_max = r_star;

      h_cells = round((y_max-y_min) / ((z_max-z_min)/v_cells));
      dh = (y_max-y_min)/h_cells;
      dv = (z_max-z_min)/v_cells;
      //cout << "dh " << dh << endl;
      //cout << "dv " << dv << endl;

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

      
      build_grid(h_grid, v_grid, h_cells, v_cells);
      grid_cells(h_grid, v_grid, patches, h_cells, v_cells);
          particles_calc.clear();
         for (int w=0; w<h_cells; w++){
         taus.push_back({});
         for (int y=0; y<v_cells; y++){
            taus[w].push_back(0.0);
         }
         }
         
         double forward_flux = 0.;
         
         for (dust& p : particles) {
              if (p.timestamp == timestamps[it]) {
               //if (p.timestamp == 0.0) {
              
                   particles_calc.push_back(p);
                 if (p.x_dp > 0.0) {
                 forward_flux = forward_flux + forward_scat(opac.stellar_gsc(p.size), 
                 opac.stellar_scat(p.size), p.m, p.x_dp*a_p, p.y_dp*a_p, p.z_dp*a_p, p.size, p.tau); }
           
              }
         }
          
         extinction(particles_calc, patches,  h_grid, v_grid,taus, h_cells, v_cells, 
         z_max, z_min, y_max, y_min);
        
         output.write((char*) &timestamps[it], sizeof(double));
         
         f_test = flux(taus, patches, h_cells, v_cells);
         cout << "extinction " << f_test << endl;
         cout << "scattering " << forward_flux << endl;
         output.write((char*) &f_test, sizeof(double));
         output.write((char*) &forward_flux, sizeof(double));
         f_test = f_test + forward_flux;
         cout << "total " << f_test << endl;
        //cout << forward_flux << endl;
         
         output.write((char*) &f_test, sizeof(double));

         if (f_test<f_test_o) {
          f_test_o = f_test;
          } 

          taus.clear();
          h_grid.clear();
          v_grid.clear();
          patches.clear();

           
      
    }

    transit_depths.push_back(f_test_o);
    cout << "transit depth " << f_test_o << endl;
    grid_cell_size.push_back(dv);


    }

   return 0;
}