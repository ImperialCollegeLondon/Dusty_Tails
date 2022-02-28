#include <iostream>
#include<iostream>
#include<fstream>
#include<vector>
#include <cmath>

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

const int x_cells = 100;
const int y_cells = 100;

double r_star = 0.24; //stellar radius in terms of the semimajor axis

double dx = (2.0*0.24)/x_cells;
double dy = (2.0*0.24)/y_cells;

double x_grid[x_cells];
double y_grid[y_cells];

double patches[x_cells][y_cells];

const int timesteps = 200;

double angle[timesteps];

double* angles(double times[timesteps], double t0, double angle[timesteps]){
   for (int i=0; i<=timesteps; i++){
      angle[i] = 2.0*M_PI * (times[i] - t0);
   }
   return angle;
}

void build_grid(double x_grid[x_cells], double y_grid[y_cells]){
   x_grid[0] = -1.0 *r_star;

   for (unsigned int i=1;i<=x_cells; i++){
      x_grid[i] = x_grid[i-1] + dx;
      
   }
   y_grid[0] = -1.0 *r_star;
   for (unsigned int k=1;k<=y_cells; k++){
      y_grid[k] = y_grid[k-1] + dy;
      
      cout << y_grid[k] << endl;
   }
}

void grid_cells(double x_grid[x_cells], double y_grid[y_cells], double patches[x_cells][y_cells]){
   double p1[2], p2[2], p3[2], p4[2];
   double r1, r2, r3, r4;
   double check;
   double delta_x, delta_y;
   double h, a, b;
   for (int m=0; m<x_cells; m++){
      for (int n=0; n<y_cells; n++){
            //cout <<  m << " " << n << " " << endl;
            check = 0;
            p1[0] = x_grid[m];
            //cout << "x1 " << p1[0] << " y1 " << p1[1] << endl;
            p1[1] = y_grid[n];
            p2[0] = x_grid[m+1];
            p2[1] = y_grid[n];
            p3[0] = x_grid[m+1];
            p3[1] = y_grid[n+1];
            p4[0] = x_grid[m];
            p4[1] = y_grid[n+1];
           
           
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
               patches[m][n] = dx*dy;
               
            }

            else if (r1 <r_star) {
               //Triangle shape, 1st quadrant, 1 in, 3 out
               if ((r2 >= r_star) && (r3>=r_star) && (r4 >=r_star) ) {
                  delta_x = pow(pow(r_star,2)-pow(p1[1],2), 0.5) - abs(p1[0]);
                  delta_y = pow(pow(r_star,2)-pow(p1[0],2), 0.5) - abs(p1[1]);
                  patches[m][n] = (delta_x * delta_y)/2.;
                  cout << "t 1st" << endl;
                  
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
   build_grid(x_grid, y_grid);
   grid_cells(x_grid, y_grid, patches);
   return 0;
}