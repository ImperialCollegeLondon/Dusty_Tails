//RK Dormand-Prince 5(4)

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "constants.h"
#incude "butcher.h"

using namespace std;

//Variable declaration

const double Mstar_kg = 0.8*Msun; //mass of star in kg
const double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
const double Period_days = 15.0/24.0; //period of dust grain in days
const double T = 15*3600.0; //period of dust grain in seconds

const double atol = 1e-6;
const double rtol = 1e-6;


double k1_xdot, k2_xdot, k3_xdot, k4_xdot, k5_xdot, k6_xdot, k7_xdot;
double k1_ydot, k2_ydot, k3_ydot, k4_ydot, k5_ydot, k6_ydot, k7_ydot;
double k1_zdot, k2_zdot, k3_zdot, k4_zdot, k5_zdot, k6_zdot, k7_zdot;

double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x;
double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y;
double k1_z, k2_z, k3_z, k4_z, k5_z, k6_z, k7_z;

double G_dim, a;


//function to evaluate ODE for velocity
double semimajor(double period) { //function to evaluate semi-major axis value, period in days, a in AU
    
   semi = pow(7.496e-6 * Mstar_sun * pow(Period_days, 2.0), (1.0/3.0)); //in AU
   return semi*AU_to_m; //output in meters
}


double vel_ODE(double G_dim, double pos, double x, double y, double z){
     
  vel_dot = (-G_dim * pos) / pow( ( pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) ) , (3.0/2.0) );
  return vel_dot;
}


//function to evaluate k values for velocity solver
void k_velocities(double h, double x, double y, double z, bool order7){
 
 k1_xdot = h*vel_ODE(G_dim, x, x, y, z);
 k1_ydot = h*vel_ODE(G_dim, y, x, y, z);
 k1_zdot = h*vel_ODE(G_dim, z, x, y, z);

 k2_xdot = h*vel_ODE(G_dim, x + a21*k1_xdot , x + a21*k1_xdot, y + a21*k1_ydot, z + a21*k1_zdot);
 k2_ydot = h*vel_ODE(G_dim, y + a21*k1_ydot , x + a21*k1_xdot, y + a21*k1_ydot, z + a21*k1_zdot);
 k2_zdot = h*vel_ODE(G_dim, z + a21*k1_zdot , x + a21*k1_xdot, y + a21*k1_ydot, z + a21*k1_zdot);

 k3_xdot = h*vel_ODE(G_dim, x + a31*k1_xdot + a32*k2_xdot , \
                        x + a31*k1_xdot + a32*k2_xdot, \
                        y + a31*k1_ydot + a32*k2_ydot , \
                        z + a31*k1_zdot + a32*k2_zdot);

 k3_ydot = h*vel_ODE(G_dim, y + a31*k1_ydot + a32*k2_ydot , \
                        x + a31*k1_xdot + a32*k2_xdot, \
                        y + a31*k1_ydot + a32*k2_ydot , \
                        z + a31*k1_zdot + a32*k2_zdot);

 k3_zdot = h*vel_ODE(G_dim, z + a31*k1_zdot + a32*k2_zdot , \
                        x + a31*k1_xdot + a32*k2_xdot, \
                        y + a31*k1_ydot + a32*k2_ydot , \
                        z + a31*k1_zdot + a32*k2_zdot);

 k4_xdot = h*vel_ODE(G_dim, x + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
                        x + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
                        y + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
                        z + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

 k4_ydot = h*vel_ODE(G_dim, y + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
                        x + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
                        y + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
                        z + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

 k4_zdot = h*vel_ODE(G_dim, z + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot, \
                        x + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
                        y + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
                        z + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);
 
 k5_xdot = h*vel_ODE(G_dim, x + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
                        x + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
                        y + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
                        z + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

 k5_ydot = h*vel_ODE(G_dim, y + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
                        x + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
                        y + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
                        z + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

 k5_zdot = h*vel_ODE(G_dim, z + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot, \
                        x + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
                        y + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
                        z + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

 k6_xdot = h*vel_ODE(G_dim, x + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
                        x + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
                        y + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
                        z + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);

 k6_ydot = h*vel_ODE(G_dim, y + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
                        x + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
                        y + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
                        z + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);

 k6_zdot = h*vel_ODE(G_dim, z + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot, \
                        x + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
                        y + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
                        z + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);

 if (order7 = TRUE) {
   
  k7_xdot = h*vel_ODE(G_dim, x + a71*k1_xdot + a72*k2_xdot + a73*k3_xdot + a74*k4_xdot + a75*k5_xdot + a76*k6_xdot, \
                         x + a71*k1_xdot + a72*k2_xdot + a73*k3_xdot + a74*k4_xdot + a75*k5_xdot + a76*k6_xdot, \
                         y + a71*k1_ydot + a72*k2_ydot + a73*k3_ydot + a74*k4_ydot + a75*k5_ydot + a76*k6_ydot, \
                         z + a71*k1_zdot + a72*k2_zdot + a73*k3_zdot + a74*k4_zdot + a75*k5_zdot + a76*k6_zdot);

  k7_ydot = h*vel_ODE(G_dim, y + a71*k1_ydot + a72*k2_ydot + a73*k3_ydot + a74*k4_ydot + a75*k5_ydot + a76*k6_ydot, \
                         x + a71*k1_xdot + a72*k2_xdot + a73*k3_xdot + a74*k4_xdot + a75*k5_xdot + a76*k6_xdot, \
                         y + a71*k1_ydot + a72*k2_ydot + a73*k3_ydot + a74*k4_ydot + a75*k5_ydot + a76*k6_ydot, \
                         z + a71*k1_zdot + a72*k2_zdot + a73*k3_zdot + a74*k4_zdot + a75*k5_zdot + a76*k6_zdot);

  k7_zdot = h*vel_ODE(G_dim, z + a71*k1_zdot + a72*k2_zdot + a73*k3_zdot + a74*k4_zdot + a75*k5_zdot + a76*k6_zdot;  \
                        x + a71*k1_xdot + a72*k2_xdot + a73*k3_xdot + a74*k4_xdot + a75*k5_xdot + a76*k6_xdot, \
                        y + a71*k1_ydot + a72*k2_ydot + a73*k3_ydot + a74*k4_ydot + a75*k5_ydot + a76*k6_ydot, \
                        z + a71*k1_zdot + a72*k2_zdot + a73*k3_zdot + a74*k4_zdot + a75*k5_zdot + a76*k6_zdot); 

   }


}

void k_positions(double h, double vx, double vy, double vz, bool order7) {
    
    k1_x = h*vx;
    k1_y = h*vy;
    k1_z = h*vz;
    
    k2_x = h* (vx + a21*k1_x);
    k2_y = h* (vy + a21*k1_y);
    k2_z = h* (vz + a21*k1_z);
    
    k3_x = h* (vx + a31*k1_x + a32*k2_x);
    k3_y = h* (vy + a31*k1_y + a32*k2_y);
    k3_z = h* (vz + a31*k1_z + a32*k2_z);
    
    k4_x = h* (vx + a41*k1_x + a42*k2_x + a43*k3_x);
    k4_y = h* (vy + a41*k1_y + a42*k2_y + a43*k3_y);
    k4_z = h* (vz + a41*k1_z + a42*k2_z + a43*k3_z);
    
    k5_x = h* (vx + a51*k1_x + a52*k2_x + a53*k3_x + a54*k4_x);
    k5_y = h* (vy + a51*k1_y + a52*k2_y + a53*k3_y + a54*k4_y);
    k5_z = h* (vz + a51*k1_z + a52*k2_z + a53*k3_z + a54*k4_z);
    
    k6_x = h* (vx + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x);
    k6_y = h* (vy + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y);
    k6_z = h* (vz + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z);
    
    if (order7 = TRUE) {
   
       k7_x = h* (vx + a71*k1_x + a72*k2_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x);
       k7_y = h* (vy + a71*k1_y + a72*k2_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y);
       k7_z = h* (vz + a71*k1_z + a72*k2_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z); 
       

   }
}

void new_velocity(double h, double vx, double vy, double vz, double x, double y, double z, bool order7){
    
    if (order7 = FALSE) {
        
        k_velocities(h, x, y, z, FALSE);
        xdot_new = vx + b1*k1_xdot + b2* k2_xdot + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;        
        ydot_new = vy + b1*k1_ydot + b2* k2_ydot + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
        zdot_new = vz + b1*k1_zdot + b2* k2_zdot + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;
        
    }
    
    else {
      k_velocities(h, x, y, z, TRUE);
      xdot_snew = vx + \
          bs1*k1_xdot + bs2* k2_xdot + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7* k7_xdot;
      ydot_snew = vy + \
          bs1*k1_ydot + bs2* k2_ydot + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7* k7_ydot;
      zdot_snew = vz + \
          bs1*k1_zdot + bs2* k2_zdot + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7* k7_zdot; 
    }
    
}


void new_position(double h, double vx, double vy, double vz, double x, double y, double z, bool order7){
    
    if (order7 = FALSE) {
     k_positions(h, vx, vy, vx, FALSE);
     xdot_new = vx + b1*k1_xdot + b2* k2_xdot + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;
     ydot_new = vy + b1*k1_ydot + b2* k2_ydot + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
     zdot_new = vz + b1*k1_zdot + b2* k2_zdot + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;    
    }
    
    else {
      k_positions(h, vx, vy, vz, TRUE);
      xdot_snew = vx + \
          bs1*k1_xdot + bs2* k2_xdot + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7* k7_xdot;
      ydot_snew = vy + \
          bs1*k1_ydot + bs2* k2_ydot + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7* k7_ydot;
      zdot_snew = vz + \
          bs1*k1_zdot + bs2* k2_zdot + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7* k7_zdot; 
    }
}

double error(double value1, double value2){
    
    scale1 = atol + abs(value1)*rtol;
    scale2 = atol + abs(value2)*rtol;
    
    delta = abs(value1 - value2);
    
    err = 
    
}



void RK_solver(double h0, vector < vector <double> > pos, vector < vector <double> > vel){
    
    for (unsigned int i = 0; i < 1000; i++) {
        
        if (i = 0) {
            
            new_velocity(h0, vel[0][0], vel[0][1], vel[0][2], pos[0][0], pos[0][1], pos[0][2], FALSE);
            new_velocity(h0, vel[0][0], vel[0][1], vel[0][2], pos[0][0], pos[0][1], pos[0][2], TRUE);
            
            new_position(h0, vel[0][0], vel[0][1], vel[0][2], pos[0][0], pos[0][1], pos[0][2], FALSE);
            new_position(h0, vel[0][0], vel[0][1], vel[0][2], pos[0][0], pos[0][1], pos[0][2], TRUE);
            
            
            
            
        }
        
        else{
            
        
        
        
            }
    }
}

int main() {
    
    a = semimajor(Period_days);
    G_dim = (G* pow(T, 2.0) * M_star_kg) / pow(a, 3.0); //dimensionless gravitational constant
    
    
}
