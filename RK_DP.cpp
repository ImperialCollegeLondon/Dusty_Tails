//RK Dormand-Prince 5(4)

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"

using namespace std;


const double Mstar_kg = 0.8*Msun; //mass of star in kg
const double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
const double Period_days = 0.8; //period of dust grain in days
const double T = 0.8*24.0*60.0*60.0; //period of dust grain in seconds

const double S = 0.9; //safety factor
const double tol = 1e-7; //error tolerance


double k1_xdot, k2_xdot, k3_xdot, k4_xdot, k5_xdot, k6_xdot, k7_xdot;
double k1_ydot, k2_ydot, k3_ydot, k4_ydot, k5_ydot, k6_ydot, k7_ydot;
double k1_zdot, k2_zdot, k3_zdot, k4_zdot, k5_zdot, k6_zdot, k7_zdot;

double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x;
double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y;
double k1_z, k2_z, k3_z, k4_z, k5_z, k6_z, k7_z;

double xdot_new, ydot_new, zdot_new;
double x_new, y_new, z_new;
double xdot_snew, ydot_snew, zdot_snew;
double x_snew, y_snew, z_snew;

double xdot4, ydot4, zdot4, x4, y4, z4;
double xdot5, ydot5, zdot5, x5, y5, z5;

double del, h_opt, max_delta, scalar;

double xdot_err, ydot_err, zdot_err;
double x_err, y_err, z_err;

double G_dim, a, h_new, h0, semi;
double vel_dot;

vector <double> deltas, delta_values, x_positions, y_positions, z_positions, time_plot, semis, t, timing, as;
vector <double> position_new, velocity_new, r0, v0, x, y, z;

vector <vector <double> > positions, velocities, final_positions, final_velocities;

double eval_velocity(double h, double vx, double vy, double vz, double x, double y, double z, bool order7);



void file_creator(vector<double> x, vector<double> y, vector<double> z);

void file_creator(vector<double> x, vector<double> y, vector<double> z) {
	ofstream myfile("test_data.txt");
	if (myfile.is_open()) {
		for (int i = 0; i < x.size(); i++) {
			char string[20];
            
			myfile << x[i] << ",";
			myfile << y[i] << ",";
            myfile << z[i] << "\n";
		}
		myfile.close();
	}
    
	else cout << "unable to open file";

}

void file_plot(vector<double> a, vector<double> b);

void file_plot(vector<double> a, vector<double> b) {
	ofstream myfile("oi_data.txt");
	if (myfile.is_open()) {
		for (int i = 0; i < b.size(); i++) {
			char string[20];
            
			myfile << a[i] << ",";
			myfile << b[i] << "\n";
		}
		myfile.close();
	}
    
	else cout << "unable to open file";

}


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
void k_velocities(double h, double vx, double vy, double vz, double x, double y, double z, bool order){
 
 k1_xdot = h*vel_ODE(G_dim, x, x, y, z);
 k1_ydot = h*vel_ODE(G_dim, y, x, y, z);
 k1_zdot = h*vel_ODE(G_dim, z, x, y, z);
    
 eval_position(c2*h, vx + a21*k1_xdot, vy + a21*k1_ydot, vz + a21*k1_zdot, x, y, z, order);

 k2_xdot = h*vel_ODE(G_dim, x_new , x_new, y_new, z_new);
 k2_ydot = h*vel_ODE(G_dim, y_new,  x_new, y_new, z_new);
 k2_zdot = h*vel_ODE(G_dim, z_new,  x_new, y_new, z_new);
    
 eval_position(c3*h, vx + a31*k1_xdot + a32*k2_xdot, \
              vy + a31*k1_ydot + a32*k2_ydot,\
              vz + a31*k1_zdot + a32*k2_zdot, order);

 k3_xdot = h*vel_ODE(G_dim, x_new, x_new, y_new, z_new);
 k3_ydot = h*vel_ODE(G_dim, y_new, x_new, y_new, z_new);
 k3_zdot = h*vel_ODE(G_dim, z_new, x_new, y_new, z_new);
    
 eval_position(c4*h, vx + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
              vy + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot,\
              vz + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot, order);

 k4_xdot = h*vel_ODE(G_dim, x_new, x_new, y_new, z_new);
 k4_ydot = h*vel_ODE(G_dim, y_new, x_new, y_new, z_new);
 k4_zdot = h*vel_ODE(G_dim, z_new, x_new, y_new, z_new);
	
 eval_position(c5*h, vx + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
	      vy + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
	      vz + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot, order); 
 
 k5_xdot = h*vel_ODE(G_dim, x_new, x_new, y_new, z_new);
 k5_ydot = h*vel_ODE(G_dim, y_new, x_new, y_new, z_new);
 k5_zdot = h*vel_ODE(G_dim, z_new, x_new, y_new, z_new);
 
 eval_position(c6*h, vx + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
	       vy + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
	       vz + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot, order);

 k6_xdot = h*vel_ODE(G_dim, x_new, x_new, y_new, z_new);
 k6_ydot = h*vel_ODE(G_dim, y_new, x_new, y_new, z_new);
 k6_zdot = h*vel_ODE(G_dim, z_new, x_new, y_new, z_new);	

 if (order == true) {
	 
  eval_position(c7*h, vx + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75*k5_xdot + a76*k6_xdot, \
		vy + a71*k1_ydot  + a73*k3_ydot + a74*k4_ydot + a75*k5_ydot + a76*k6_ydot, \
		vz + a71*k1_zdot  + a73*k3_zdot + a74*k4_zdot + a75*k5_zdot + a76*k6_zdot, order);
   
  k7_xdot = h*vel_ODE(G_dim, x_new, x_new, y_new, z_new);
  k7_ydot = h*vel_ODE(G_dim, y_new, x_new, y_new, z_new);
  k7_zdot = h*vel_ODE(G_dim, z_new, x_new, y_new, z_new); 

 }


}


void k_positions(double h, double vx, double vy, double vz, double x, double y, double z, bool order) {
    
    k1_x = h*vx;
    k1_y = h*vy;
    k1_z = h*vz;
    
    eval_velocity(c2*h, vx, vy, vz, x + a21*k1_x, y + a21*k1_y, z + a21*k1_z, order);
    
    k2_x = h* xdot_new;
    k2_y = h* ydot_new;
    k2_z = h* zdot_new;
    
    eval_velocity(c3*h, vx, vy, vz, \
             x + a31*k1_x + a32*k2_x, \
             y + a31*k1_y + a32*k2_y, \
             z + a31*k1_z + a32*k2_z , order);
    
    k3_x = h* xdot_new;
    k3_y = h* ydot_new;
    k3_z = h* zdot_new;
    
    eval_velocity(c4*h, vx, vy, vz, \
             x + a41*k1_x + a42*k2_x + a43*k3_x, \
             y + a41*k1_y + a42*k2_y + a43*k3_y, \
             z + a41*k1_z + a42*k2_z +  a43*k3_z , order);
    
    
    k4_x = h* xdot_new;
    k4_y = h* ydot_new;
    k4_z = h* zdot_new;
    
    eval_velocity(c5*h, vx, vy, vz, \
             x + a51*k1_x + a52*k2_x + a53*k3_x + a54*k4_x, \
             y + a51*k1_y + a52*k2_y + a53*k3_y + a54*k4_y, \
             z + a51*k1_z + a52*k2_z + a53*k3_z + a54*k4_z , order);
    
    k5_x = h* xdot_new;
    k5_y = h* ydot_new;
    k5_z = h* zdot_new;
    
    eval_velocity(c6*h, vx, vy, vz, \
             x + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
             y + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
             z + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z , order);
    
    k6_x = h* xdot_new;
    k6_y = h* ydot_new;
    k6_z = h* zdot_new;
    
    
    if (order == true) {
        
       eval_velocity(c7*h, vx, vy, vz, \
             x + a71*k1_x  + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
             y + a71*k1_y  + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
             z + a71*k1_z  + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z , order);
   
       k7_x = h* xdot_new;
       k7_y = h* ydot_new;
       k7_z = h* zdot_new;
       
    }
}



double eval_velocity(double h, double vx, double vy, double vz, double x, double y, double z, bool order){
    
    if (order == false) {
        
        k_velocities(h, x, y, z, false);
        
        xdot_new = vx + b1*k1_xdot  + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;
        ydot_new = vy + b1*k1_ydot  + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
        zdot_new = vz + b1*k1_zdot  + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;
        
    } else {
        
       k_velocities(h, x, y, z, true);
        
       xdot_new = vx + \
          bs1*k1_xdot  + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7* k7_xdot;
       ydot_new = vy + \
          bs1*k1_ydot  + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7* k7_ydot;
       zdot_new = vz + \
          bs1*k1_zdot  + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7* k7_zdot; 
    }
    
    return xdot_new, ydot_new, zdot_new;
}


double eval_position(double h, double vx, double vy, double vz, double x, double y, double z, bool order){
    
    if (order == false) {
     k_positions(h, vx, vy, vx, x, y, z, false);
        
     x_new = x + b1*k1_x + b3* k3_x + b4* k4_x + b5* k5_x + b6* k6_x;
     y_new = y + b1*k1_y + b3* k3_y + b4* k4_y + b5* k5_y + b6* k6_y;
     z_new = z + b1*k1_z + b3* k3_z + b4* k4_z + b5* k5_z + b6* k6_z;    
    
    } else {
        
      k_positions(h, vx, vy, vz, x, y, z, true);
        
      x_new = x + bs1*k1_x + bs3* k3_x + bs4* k4_x + bs5* k5_x + bs6* k6_x + bs7 * k7_x;
      y_new = y + bs1*k1_y + bs3* k3_y + bs4* k4_y + bs5* k5_y + bs6* k6_y + bs7 * k7_y;
      z_new = z + bs1*k1_z + bs3* k3_z + bs4* k4_z + bs5* k5_z + bs6* k6_z + bs7 * k7_z;
    }
    
    return x_new, y_new, z_new;
}



double delta( double value1, double value2){
        
    del = fabs(value1 - value2) + 1e-10;
   
    
    return del/(tol);
    
      
}

double h_optimal(vector <double> deltas, double h){
    
    max_delta = *max_element(deltas.begin(), deltas.end()); //includes tolerance
    
    h_opt = S * h * pow(1/max_delta, 1.0/5.0);
    
    if (h_opt < h){
        
        return h_opt;
        
    } else return h;

}




vector <double> h_check(double h, vector <double> velocity, vector <double> position){
    
    deltas.clear();
    
    //4th order velocity
    eval_velocity(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], false);
    //cout << "4th order " << xdot_new << endl;
    xdot4 = xdot_new;
    ydot4 = ydot_new;
    zdot4 = zdot_new;

    //5th order velocity
    eval_velocity(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], true);
    //cout << "5th order " << xdot_new << endl;
    xdot5 = xdot_new;
    ydot5 = ydot_new;
    zdot5 = zdot_new;
    
    //4th order position
    eval_position(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], false);
    
    x4 = x_new;
    y4 = y_new;
    z4 = z_new; 
    
    //5th order postiion
    eval_position(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], true);
    
    x5 = x_new;
    y5 = y_new;
    z5 = z_new;
    
    //error on xdot
    xdot_err = delta(xdot4, xdot5);
    deltas.push_back(xdot_err);
    
    //error on ydot
    ydot_err = delta(ydot4, ydot5);
    deltas.push_back(ydot_err);
    
    //error on zdot
    zdot_err = delta(zdot4, zdot5);
    deltas.push_back(zdot_err);
    
    
    //error on x  
    x_err = delta(x4, x5);
    deltas.push_back(x_err);
    
    //error on y
    y_err = delta(y4, y5);
    deltas.push_back(y_err);
    
    //error on z
    z_err = delta(z4, z5);
    deltas.push_back(z_err);
    
    return deltas;
 
}


vector <double> new_variables(double h, vector <double> velocity, vector <double> position) {
    
    //ensure vector are clear
    velocity_new.clear();
    position_new.clear();
    
    eval_velocity(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], false); 
    
    velocity_new.push_back(xdot_new);
    velocity_new.push_back(ydot_new);
    velocity_new.push_back(zdot_new);
       
        
    eval_position(h, velocity[0], velocity[1], velocity[2], position[0], position[1], position[2], false);
    
     
    position_new.push_back(x_new);
    position_new.push_back(y_new);
    position_new.push_back(z_new);
         
    return velocity_new;
    return position_new;
}





void RK_solver(double h0, vector < vector <double> > pos, vector < vector <double> > vel, \
     vector <double> x_positions, vector <double> y_positions, vector <double> z_positions, vector <double> t){
    
    
    //obtain delta values for the 6 variables
    delta_values = h_check(h0, vel[0], pos[0]);
        
        
        //calculate new h
    h_new = h_optimal(delta_values, h0);
        
    new_variables(h_new, vel[0], pos[0]);
        
    x_positions.push_back(position_new[0]);
    y_positions.push_back(position_new[1]);
    z_positions.push_back(position_new[2]);
        
    pos.push_back(position_new);
    vel.push_back(velocity_new);
        
        
    t.push_back(h_new);
        
    
    for (unsigned int i = 1; i < 1e+6; i++) {
        
        //obtain delta values for the 6 variables
        delta_values = h_check(h_new, vel[i], pos[i]);
        
        
        //calculate new h
        h_new = h_optimal(delta_values, h_new);
        t.push_back(t[i] + h_new);
        //cout << h_new << endl;
        
        new_variables(h_new, vel[i], pos[i]);
        
        x_positions.push_back(position_new[0]);
        y_positions.push_back(position_new[1]);
        z_positions.push_back(position_new[2]);
        
        pos.push_back(position_new);
        vel.push_back(velocity_new);
        
        
        
       
       
     } 
    
    
     final_positions = pos;
     final_velocities = vel;
     x = x_positions;
     y = y_positions;
     z = z_positions;
     
     timing = t;
     
     cout << "done!" << endl;
}

vector <double> semimajor_plot(vector <double> x_pos, vector <double> y_pos, vector <double> z_pos){
     
     for (unsigned int i = 0; i < 1e+6; i++) {
     
         scalar = pow( pow(x_pos[i], 2.) + pow(y_pos[i], 2.) + pow(z_pos[i], 2.), 0.5 );
     
         semis.push_back(scalar);
     
     }
     
     return semis;
}

int main() {
    
    
    a = semimajor(Period_days);
    G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant
    
    h0 = 0.001;
    
    //Define initial position in dimensionless units
    
    double x0 = 1.0;
    double y0 = 0.0;
    double z0 = 0.0;
    
    //initial position vector
    r0.push_back(x0);
    r0.push_back(y0);
    r0.push_back(z0);
    
    //Define initial velocity in dimensionless units
    
    double xdot0 = 0.0;
    double ydot0 = 2.0*PI;
    double zdot0 = 0.0; 
    
    //initial velocity vector
    v0.push_back(xdot0);
    v0.push_back(ydot0);
    v0.push_back(zdot0);
    
    positions.push_back(r0);
    velocities.push_back(v0);
    
    x_positions.push_back(x0);
    y_positions.push_back(y0);
    z_positions.push_back(z0);
    
    double t0 = 0.0;
    
    t.push_back(t0);

    RK_solver(h0, positions, velocities, x_positions, y_positions, z_positions, t);
    
    as = semimajor_plot(x, y, z);
    
    file_plot(timing, as);
    
    file_creator(x, y, z);
       
    
}
