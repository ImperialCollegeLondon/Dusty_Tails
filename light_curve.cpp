#include <iostream>
#include<iostream>
#include<fstream>
#include<vector>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

#define amu 1.661e-24
#define kb 1.381e-16

#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define Msun_cgs 1.9885e+33 // solar mass
#define Mearth_cgs 5.972e+27 // earth mass grams

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
};


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
double mdot_read = 1.0; //mass loss rate of planet in Earth masses per Gyr
double mdot = mdot_read*Mearth_cgs/gyr;
double s_0 = 1e-4;
double rho_d = 4.0;
double T = 15.68*60.*60.; //planetary period in seconds
double mbig = (mdot * T * 0.01) / 500. ; // 0.01 dependent on when particles are being thrown out of planet
double n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
double m_star = 0.66; //stellar mass in solar masses
double a_p = pow((G_cgs*m_star*Msun_cgs* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in cgs
double inclination = 1.425; //in radians
int h_cells, v_cells;

double r_star = 0.2415; //stellar radius in terms of the semimajor axis

double z_min, z_max, y_min, y_max;
double dh, dv; 

const int timesteps = 100;

double angle[timesteps];

ofstream output("./light_curves/cellsize_vs_transitdepth.bin", ios::out | ios::binary); 

double* angles(double times[timesteps], double t0, double angle[timesteps]){
   for (int i=0; i<=timesteps; i++){
      angle[i] = 2.0*M_PI * (times[i] - t0);
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
  std::fstream output;
  output.open("./simulations/KIC1255b_Al2O3_1micro_1mdot_sph_1orb/output_struct_format.bin", std::fstream::in | std::fstream::binary);
  output.seekg(0, ios::end);
  int size=output.tellg();
  output.seekg(0, ios::beg);
  //cout << "size " << size << endl;
  long int total = size/sizeof(dust);
 // cout << "total " << total << endl;
  vector <dust_read> dust_grains_out;
  for(int i = 0; i < total; i++){
       dust_grains_out.push_back(dust_read());
       output.read((char *) &dust_grains_out[i], sizeof(dust_read));
       
       if (dust_grains_out[i].id == 0) {
          //cout << "id is zero " << endl;
          //dust_grains_out.erase(dust_grains_out.end());
          dust_grains_out.pop_back();
          //cout << "going to break " << endl;
          break;
          
       }
       //cout << dust_grains_out[i].id << endl;
    }
   output.close();
   return dust_grains_out;
}

vector <double> scaled_pos(double x, double y, int h_cells, int v_cells, 
double z_max, double z_min, double y_max, double y_min){
   double x_scaled, y_scaled;
   x_scaled = (h_cells/(y_max-y_min))*x - ((h_cells*y_min)/(y_max-y_min));
   y_scaled = (v_cells/(z_max-z_min))*y - ((v_cells*z_min)/(z_max-z_min));
   //cout << "x " << x << " x_scaled " << x_scaled << endl;
   return {x_scaled, y_scaled};
}

void  extinction( vector <dust> particles, vector <vector <double>> &patches, 
         vector<double> &h_grid, vector<double> &v_grid,
         vector<vector <double>> &taus, int h_cells, int v_cells,
         double z_max, double z_min, double y_max, double y_min){
         
         //cout << "At extinction " << endl;
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
                     //cout << "h " << h_index[i] << endl;
                     //cout << "v " << v_index[j] << endl;
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
            // cout << "r1 " << r1 << endl;
            // cout << "r2 " << r2 << endl;
            // cout << "r3 " << r3 << endl;
            // cout << "r4 " << r4 << endl;
            
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

        f = f + ((patches[m][n] * exp(-1.0*taus[m][n])) / total_grid_area);
         //f = f+(patches[m][n]/ (PI * pow(r_star,2.)));
      }
   }
   //
   return f;
}

int main(){

using namespace std;
#include <cstring>
   vector <double> period_steps = {};
   vector <dust_read> particles_read;
   vector <dust> particles;
   particles_read = read_data();
   double t0 = 0.0;
   int counter = 0;
   double f_test_o = 1.0;
   vector <double> timestamps;
   vector <double> transit_depths = {};
   vector <double> grid_cell_size = {};

   for( dust_read& p : particles_read) {
       particles.push_back(dust());
      
       particles[counter].timestamp = p.timestamp;
       double key = p.timestamp;
       if (find(timestamps.begin(), timestamps.end(), key) == timestamps.end()) {
          timestamps.push_back(key);
      }

       particles[counter].phi = 2.0*M_PI * (p.timestamp - t0);
       particles[counter].x_p = p.x_dust * cos(particles[counter].phi) - p.y_dust * sin(particles[counter].phi);
       particles[counter].y_p = p.x_dust * sin(particles[counter].phi) + p.y_dust * cos(particles[counter].phi);
       particles[counter].z_p = p.z_dust;
       particles[counter].x_dp = particles[counter].x_p * sin(inclination) - particles[counter].z_p * cos(inclination);
       particles[counter].y_dp = particles[counter].y_p;
       particles[counter].z_dp = particles[counter].x_p * cos(inclination) + particles[counter].z_p * sin(inclination);
       particles[counter].m = p.m_dust;
       particles[counter].kappa = p.kappa_planck;
       
       counter = counter + 1;
       }
     
      
   vector< double > h_grid;
   vector<double> v_grid;
   vector<vector<double>> patches;
   vector<vector<double>> taus;
   vector <dust> particles_calc;
   //cout << "n_mini " << n_mini << endl;
   for (int i=0; i<1; i++){
      
     // cout << "i " << i << endl;
      //memset(taus, 0, sizeof(taus));
      v_cells = i+5;

      double f_test = 0.;
      double f_test_o = 1.0;

      for ( int m=0; m<timestamps.size(); m++) {
      double theta;
      double z_planet;
      //cout << "timestamps[m] " << endl;
      theta = 2.0*PI*timestamps[m];
      //cout << "theta " << theta <<  endl;
      z_planet = cos(theta) * cos(inclination);
      //cout << "z_planet " << z_planet << endl;
      z_min = z_planet - 0.02;
      z_max = z_planet + 0.02;
      y_min = -sqrt(pow(r_star,2)- pow(z_min,2));
      y_max = -1.0*y_min;

      h_cells = round((y_max-y_min) / ((z_max-z_min)/v_cells));
      //cout << "h_cells " << h_cells << endl;
      //cout << "v_cells " << v_cells << endl;
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
         //cout << "patches m" << patches[m].size() << endl;
      }
      //cout << "patches " << patches.size() << endl;
      
      build_grid(h_grid, v_grid, h_cells, v_cells);
      grid_cells(h_grid, v_grid, patches, h_cells, v_cells);
         particles_calc.clear();
         for (int w=0; w<h_cells; w++){
         taus.push_back({});
         for (int y=0; y<v_cells; y++){
            taus[w].push_back(0.0);
         }
         }
         //cout << "taus " << taus.size() << endl;
         for (dust& p : particles) {
            if (p.timestamp == timestamps[m]) {
               particles_calc.push_back(p);
            }
         }
         extinction(particles_calc, patches,  h_grid, v_grid,taus, h_cells, v_cells, 
         z_max, z_min, y_max, y_min);
        
         f_test = flux(taus, patches, h_cells, v_cells);
         cout << f_test << endl;
         if (f_test<f_test_o) {
          f_test_o = f_test;
         } 

          taus.clear();
          h_grid.clear();
          v_grid.clear();
          patches.clear();
      
    }

    transit_depths.push_back(f_test_o);
    cout << "depth " << f_test_o << endl;
    grid_cell_size.push_back(dv);


    }

   for (int i =0; i<50; i++){
       output.write((char*) &grid_cell_size[i], sizeof(double));
       output.write((char*) &transit_depths[i], sizeof(double));
   }
   return 0;
}