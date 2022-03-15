#include <iostream>
#include<iostream>
#include<fstream>
#include<vector>
#include <cmath>
#include <algorithm>

using namespace std;

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
   double kappa_dust;
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
    double p_size; //size of the particle
    double p_opacity; //particle efficiency coeff
    double p_mass; //mass of particle
    double p_tau; //optical depth
    double p_temp; //particle temperature

};
#define PI 3.14159
const double gyr = pow(10.,9) * 365. * 24. * 60. * 60.;
#define Mearth_cgs 5.972e+27 // earth mass grams
double mdot_read = 1.0; //mass loss rate of planet in Earth masses per Gyr
double mdot = mdot_read*Mearth_cgs/gyr;
double s_0 = 3e-4;
double rho_d = 4.0;
double T = 15.68*60.*60.; //planetary period in seconds
double mbig = (mdot * T * 0.01) / 1000. ; // 0.01 dependent on when particles are being thrown out of planet
double n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
double sm_to_cm = 1.93e+11; 
double inclination = 1.43; //in radians
const int h_cells = 100;
const int v_cells = 100;

double r_star = 0.24; //stellar radius in terms of the semimajor axis

double dh = (2.0*0.24)/h_cells;
double dv = (2.0*0.24)/v_cells;



double h_grid[h_cells+1];
double v_grid[v_cells+1];

double patches[h_cells+1][v_cells+1];
double taus[h_cells+1][v_cells+1];

const int timesteps = 200;

double angle[timesteps];

double* angles(double times[timesteps], double t0, double angle[timesteps]){
   for (int i=0; i<=timesteps; i++){
      angle[i] = 2.0*M_PI * (times[i] - t0);
   }
   return angle;
}

void build_grid(double h_grid[h_cells+1], double v_grid[v_cells+1]){
   h_grid[0] = -1.0 *r_star;
   for (unsigned int i=1;i<=h_cells; i++){
      h_grid[i] = h_grid[i-1] + dh;  
     
   }
   v_grid[0] = -1.0 *r_star;
   for (unsigned int k=1;k<=v_cells; k++){
      v_grid[k] = v_grid[k-1] + dv;
   }
}

vector <dust> read_data(){
  std::fstream output;
  output.open("./simulations/K222b_03micro_1mdot_day_3orb_tc.bin", std::fstream::in | std::fstream::binary);
  output.seekg(0, ios::end);
  int size=output.tellg();
  output.seekg(0, ios::beg);
  cout << "size " << size << endl;
  long int total = size/sizeof(dust);
  cout << "total " << total << endl;
  vector <dust> dust_grains_out;
  for(int i = 0; i < total; i++){
       dust_grains_out.push_back(dust());
       output.read((char *) &dust_grains_out[i], sizeof(dust));
    }
   output.close();
   return dust_grains_out;
}

vector <double> scaled_pos(double x, double y){
   double x_scaled, y_scaled;
   x_scaled = (h_cells/(2.0*r_star))*x + (h_cells/2.);
   y_scaled = (v_cells/(2.0*r_star))*y + (v_cells/2.);

   return {x_scaled, y_scaled};
}

void idk( vector <dust> particles, double patches[h_cells+1][v_cells+1], 
         double h_grid[h_cells+1], double v_grid[v_cells+1],
         double taus[h_cells+1][v_cells+1]){
         cout << "at idk " << endl;
         for ( dust& p : particles) {
         double rp, hp, vp, dA, dA_cell, sigma, shadow;
         vector <double> spos;
         vector <double> h_deltas, v_deltas;
         vector <int> h_index, v_index;
         
         int hit, vit;
         hp = p.y_dp;
         vp = p.z_dp;
         rp = pow(pow(hp, 2.0) + pow(vp, 2.0), 0.5);
         cout << hp << endl;
         cout << vp << endl;
         cout << rp << endl;
         cout << r_star << endl;
         if (rp < r_star) { cout << "true " << endl;}
         if (rp < r_star) {
            cout << "hp " << hp << endl;
            cout << "vp " << vp << endl;
            spos = scaled_pos(hp, vp);
            cout << "hp s " << spos[0] << endl;
            cout << "vp s " << spos[1] << endl;
            hit = floor(spos[0]);
            vit = floor(spos[1]);

            // if ((hit == 0) && (hp-(dh/2.) < h_grid[0])){
            //       h_index = {-1, 0};
            //       h_deltas[0] = 0.0;
            //       h_deltas[1] = abs(abs(hp+(dh/2.)) -abs(h_grid[0]));
            // }
            // else if ((hit == (h_cells -1)) && (hp+(dh/2.) > h_grid[h_cells])) {
            //    h_index = {h_cells, -1};
            //    h_deltas[0] = abs( abs(h_grid[h_cells]) - abs(hp-(dh/2.)));
            //    h_deltas[1] = 0.0;
            // }

            // else if (hp > (h_grid[hit]+(dh/2.)) ) {
            //    h_index = {hit, hit+1};
            //    h_deltas[0] = abs( abs(h_grid[hit+1]) - abs(hp-(dh/2.)));
            //    h_deltas[1] = abs( abs(hp+(dh/2.)) - abs(h_grid[hit+1]));

            // } 

            // else if (hp < (h_grid[hit]+(dh/2.))) {
            //    h_index = {hit-1, hit};
            //    h_deltas[0] = abs( abs(h_grid[hit]) - abs(hp - (dh/2.)));
            //    h_deltas[1] = abs( abs(hp+(dh/2.)) - abs(h_grid[hit]));

            // }
            
            // if ((vit == 0) && (vp-(dv/2.) < v_grid[0])){
            //       v_index = {-1, 0};
            //       v_deltas[0] = 0.0;
            //       v_deltas[1] = abs(abs(vp+(dv/2.)) -abs(v_grid[0]));
            // }
            // else if ((vit == (v_cells -1)) && (vp+(dv/2.) > v_grid[v_cells])) {
            //    v_index = {v_cells, -1};
            //    v_deltas[0] = abs( abs(v_grid[v_cells]) - abs(vp-(dv/2.)));
            //    v_deltas[1] = 0.0;
            // }

            // else if (vp > (v_grid[vit]+(dv/2.)) ) {
            //    v_index = {vit, vit+1};
            //    v_deltas[0] = abs( abs(v_grid[vit+1]) - abs(vp-(dv/2.)));
            //    v_deltas[1] = abs( abs(vp+(dv/2.)) - abs(v_grid[vit+1]));

            // } 

            // else if (vp < (v_grid[vit]+(dv/2.))) {
            //    v_index = {vit-1, vit};
            //    v_deltas[0] = abs( abs(v_grid[vit]) - abs(vp - (dv/2.)));
            //    v_deltas[1] = abs( abs(vp+(dv/2.)) - abs(v_grid[vit]));

            // }
            
            // for (unsigned int i=0; i<2; i++){
            //    for (unsigned int j=0; j<2; j++){
            //       if ((h_index[i]>0) && (v_index[j]>0)) {
            //          dA = h_deltas[i]*v_deltas[j];
            //          dA_cell = patches[h_index[i]][v_index[j]];
            //          sigma = (p.kappa * p.m) / (pow(sm_to_cm, 2.0));
            //          shadow = ((dA/(dh*dv)) * sigma*n_mini) / dA_cell;
            //          taus[h_index[i]][v_index[j]] = taus[h_index[i]][v_index[j]] + shadow;
            //       }

            //    }
            //       }
               }
            }
            
         }




void grid_cells(double h_grid[h_cells+1], double v_grid[v_cells+1], double patches[h_cells+1][v_cells+1]){
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
            if ((r1 >=r_star) && (r2 >= r_star) && (r3 >=r_star) && (r4 >= r_star)){
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


int main(){
   vector <double> period_steps = {};
   vector <dust_read> particles_read;
   vector <dust> particles;
   particles = read_data();
   double t0 = 0.0;
   int counter = 0;
   for( dust_read& p : particles_read) {
       particles.push_back(dust());
       //cout << p.timestamp << endl;
       //if (std::find(period_steps.begin(), period_steps.end(), p.timestamp) == period_steps.end()) {
        //cout << p.timestamp << endl;
        //period_steps.push_back(p.timestamp);
       //} 
       particles[counter].timestamp = p.timestamp;
       particles[counter].phi = 2.0*M_PI * (p.timestamp - t0);
       particles[counter].x_p = p.x_dust * cos(particles[counter].phi) - p.y_dust * sin(particles[counter].phi);
       particles[counter].y_p = p.x_dust * sin(particles[counter].phi) + p.y_dust * cos(particles[counter].phi);
       particles[counter].z_p = p.z_dust;
       particles[counter].x_dp = particles[counter].x_p * sin(inclination) - particles[counter].z_p * cos(inclination);
       particles[counter].y_dp = particles[counter].y_p;
       particles[counter].z_dp = particles[counter].x_p * cos(inclination) + particles[counter].z_p * sin(inclination);
       particles[counter].m = p.m_dust;
       particles[counter].kappa = p.kappa_dust;
       counter = counter + 1;
       }
   cout << particles[2345].y_dp << endl;
   build_grid(h_grid, v_grid);
   grid_cells(h_grid, v_grid, patches);
   idk(particles, patches,  h_grid, v_grid,taus);
   
   return 0;
}