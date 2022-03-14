#include <iostream>
#include<iostream>
#include<fstream>
#include<vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct dust {
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

const int x_cells = 100;
const int y_cells = 100;

double r_star = 0.24; //stellar radius in terms of the semimajor axis

double dx = (2.0*0.24)/x_cells;
double dy = (2.0*0.24)/y_cells;



double x_grid[x_cells+1];
double y_grid[y_cells+1];

double patches[x_cells+1][y_cells+1];

const int timesteps = 200;

double angle[timesteps];

double* angles(double times[timesteps], double t0, double angle[timesteps]){
   for (int i=0; i<=timesteps; i++){
      angle[i] = 2.0*M_PI * (times[i] - t0);
   }
   return angle;
}

void build_grid(double x_grid[x_cells+1], double y_grid[y_cells+1]){
   x_grid[0] = -1.0 *r_star;
   for (unsigned int i=1;i<=x_cells; i++){
      x_grid[i] = x_grid[i-1] + dx;  
      //cout << i <<" "<< x_grid[i] << endl;
   }
   y_grid[0] = -1.0 *r_star;
   for (unsigned int k=1;k<=y_cells; k++){
      y_grid[k] = y_grid[k-1] + dy;
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
       //cout << dust_grains_out[i].timestamp << endl;
       //cout << dust_grains_out[i].id << endl;
       //cout << dust_grains_out[i].x_dust << endl;
       //cout << dust_grains_out[i].kappa_dust << endl;
    }
   output.close();
   return dust_grains_out;
}

vector <double> scaled_pos(double x, double y){
   double x_scaled, y_scaled;
   x_scaled = (x_cells/r_star)*x + (x_cells/2.);
   y_scaled = (y_cells/r_star)*y + (y_cells/2.);

   return {x_scaled, y_scaled};
}

void idk(double xp, double yp, double patches[x_cells+1][y_cells+1], 
         double x_grid[x_cells+1], double y_grid[y_cells+1]){

         double rp;
         vector <double> spos;
         vector <double> x_deltas, y_deltas;
         vector <int> x_index, y_index;
         int xit, yit;

         rp = pow(pow(xp, 2.0) + pow(yp, 2.0), 0.5);
         if (rp < r_star) {
            spos = scaled_pos(xp, yp);
            xit = floor(spos[0]);
            yit = floor(spos[1]);

            if ((xit == 0) && (xp-(dx/2.) < x_grid[0])){
                  x_index = {-1, 0};
                  x_deltas[0] = 0.0;
                  x_deltas[1] = abs(abs(xp+(dx/2.)) -abs(x_grid[0]));
            }
            else if ((xit == (x_cells -1)) && (xp+(dx/2.) > x_grid[x_cells])) {
               x_index = {x_cells, -1};
               x_deltas[0] = abs( abs(x_grid[x_cells]) - abs(xp-(dx/2.)));
               x_deltas[1] = 0.0;
            }

            else if (xp > (x_grid[xit]+(dx/2.)) ) {
               x_index = {xit, xit+1};
               x_deltas[0] = abs( abs(x_grid[xit+1]) - abs(xp-(dx/2.)));
               x_deltas[1] = abs( abs(xp+(dx/2.)) - abs(x_grid[xit+1]));

            } 

            else if (xp < (x_grid[xit]+(dx/2.))) {
               x_index = {xit-1, xit};
               x_deltas[0] = abs( abs(x_grid[xit]) - abs(xp - (dx/2.)));
               x_deltas[1] = abs( abs(xp+(dx/2.)) - abs(x_grid[xit]));

            }
            
            if ((yit == 0) && (yp-(dy/2.) < y_grid[0])){
                  y_index = {-1, 0};
                  y_deltas[0] = 0.0;
                  y_deltas[1] = abs(abs(yp+(dy/2.)) -abs(y_grid[0]));
            }
            else if ((yit == (y_cells -1)) && (yp+(dy/2.) > y_grid[y_cells])) {
               y_index = {y_cells, -1};
               y_deltas[0] = abs( abs(y_grid[y_cells]) - abs(yp-(dy/2.)));
               y_deltas[1] = 0.0;
            }

            else if (yp > (y_grid[yit]+(dy/2.)) ) {
               y_index = {yit, yit+1};
               y_deltas[0] = abs( abs(y_grid[yit+1]) - abs(yp-(dy/2.)));
               y_deltas[1] = abs( abs(yp+(dy/2.)) - abs(y_grid[yit+1]));

            } 

            else if (yp < (y_grid[yit]+(dy/2.))) {
               y_index = {yit-1, yit};
               y_deltas[0] = abs( abs(y_grid[yit]) - abs(yp - (dy/2.)));
               y_deltas[1] = abs( abs(yp+(dy/2.)) - abs(y_grid[yit]));

            }
            
            }
         }




void grid_cells(double x_grid[x_cells+1], double y_grid[y_cells+1], double patches[x_cells+1][y_cells+1]){
   double p1[2], p2[2], p3[2], p4[2];
   double r1, r2, r3, r4;
   double check;
   double delta_x, delta_y;
   double h, a, b;
   cout <<x_grid[100] << endl;
   for (int m=0; m<x_cells; m++){
      for (int n=0; n<y_cells; n++){
            //cout <<  m << " " << n << " " << endl;
            check = 0;
            p1[0] = x_grid[m];
            //cout << "x1 " << p1[0] << " y1 " << p1[1] << endl;
            p1[1] = y_grid[n];
            p2[0] = x_grid[m+1];

            p2[1] = y_grid[n];
             //cout << "x2 " << p2[0] << " y2 " << p2[1] << endl;
            p3[0] = x_grid[m+1];
            p3[1] = y_grid[n+1];
             //cout << "x3 " << p3[0] << " y3 " << p3[1] << endl;
            p4[0] = x_grid[m];
            p4[1] = y_grid[n+1];
             //cout << "x4 " << p4[0] << " y4 " << p4[1] << endl;
           
           
            r1 = pow(pow(p1[0],2) + pow(p1[1],2), 0.5);
            r2 = pow(pow(p2[0],2) + pow(p2[1],2), 0.5);
            r3 = pow(pow(p3[0],2) + pow(p3[1],2), 0.5);
            r4 = pow(pow(p4[0],2) + pow(p4[1],2), 0.5);

            
            //cell outside
            if ((r1 >=r_star) && (r2 >= r_star) && (r3 >=r_star) && (r4 >= r_star)){
               patches[m][n] = 0.0;
               // cout << "cell outside " << endl;
               // cout << "r1 " << r1 << endl ;
               // cout << "r2 " << r2 << endl;
               // cout << "r3 " << r3 << endl;
               // cout << "r4 " << r4 << endl;
            }

            //cell inside
            else if ((r1 < r_star) && (r2 < r_star) && (r3 < r_star) && (r4 < r_star)) {
               patches[m][n] = dx*dy;
               
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
                  patches[m][n] = (dx*dy) - 0.5*delta_x * delta_y;
                
               }
            
               //2nd quadrant, 3 in, 1 out
               else if ((r2 < r_star) && (r3 < r_star) && (r4 >=r_star)) {
                  delta_x = abs(p4[0]) - pow(pow(r_star,2)-pow(p4[1],2), 0.5);
                  delta_y = abs(p4[1]) - pow(pow(r_star,2)-pow(p4[0],2), 0.5);
                  patches[m][n] = (dx*dy) - 0.5*delta_x * delta_y;
                  
               }
               //4th quadrant, 3 in, 1 out
               else if ((r3 < r_star) && (r4 < r_star) && (r2 >= r_star)) {
                  delta_x = abs(p2[0]) - pow(pow(r_star,2)-pow(p2[1],2), 0.5);
                  delta_y = abs(p2[1]) - pow(pow(r_star,2)-pow(p2[0],2), 0.5);
                  patches[m][n] = (dx*dy) - 0.5*delta_x * delta_y;
                  
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
                  patches[m][n] = (dx*dy) - 0.5*delta_x * delta_y;
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
   vector <dust> particles;
   particles = read_data();
    for( dust& p : particles) {
       //cout << p.timestamp << endl;
       if (std::find(period_steps.begin(), period_steps.end(), p.timestamp) == period_steps.end()) {
        cout << p.timestamp << endl;
        period_steps.push_back(p.timestamp);
       } 
       }
       
   build_grid(x_grid, y_grid);
   grid_cells(x_grid, y_grid, patches);
   cout << dx << endl;
   cout << patches[0][50] << endl;
   return 0;
}