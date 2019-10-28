//First test with Leap-Frog integration on a simple analytical scenario
//Scenario is star with orbiting dust grain

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "constants.h"

using namespace std;
 
//Declare variables 
const double Mstar_kg = 0.8*Msun; //mass of star in kg
const double Mstar_sun = 0.8; //mass of star in terms of mass of the sun
const double Period_days = 0.8; //period of dust grain in days
const double T = 0.8*24.0*60.0*60.0; //period of dust grain in seconds

const double del_t = 1.0/100.0; //time step of integration


double a, grav, semi; //semi-major axis, period
double gravity; //gravity scalar
double accel_new, r_new, vel_new; //quantities after solver

double gravity_0; 
double x_new, y_new, z_new;
double vxdot_next, vydot_next, vzdot_next;
double xdot_new, ydot_new, zdot_new;

vector <double> centrifugal, c_P; //centrifugal force pre term
vector <double> r0, v0, a0, omega_0; //initial position, velocity and acceleration vectors

vector <double> X, Y, Z;
vector <double> x_positions, y_positions, z_positions;
vector <double> X_positions, Y_positions, Z_positions;

vector <double> position_new, velocity_new, acceleration_new;

vector <vector <double> > positions, velocities, accelerations;
vector <vector <double> > pos, vel, accel;
vector <vector <double> > POSITIONS, VELS, ACCELS;

void file_creator(vector<double> x, vector<double> y, vector<double> z);

//end of declarations of variables


//functions
double semimajor(double period) { //function to evaluate semi-major axis value, period in days, a in AU
    
   semi = pow(7.496e-6 * Mstar_sun * pow(Period_days, 2.0), (1.0/3.0)); //in AU
   return semi*AU_to_m; //output in meters
}

double grav_scalar(double x, double y, double z) { //gravity scalar - constant in the case of the circle, good test
    gravity =  (pow(T, 2.0) * Mstar_kg * G) / ((pow(a, 3.0) )*(pow( pow(x,2.0) + pow(y, 2.0) + pow(z, 2.0), 3.0/2.0 )));
    return gravity;
}

vector <double> cross_product(vector <double> vector_1, vector <double> vector_2) {
    
    c_P[0] = vector_1[1]*vector_2[2] - vector_1[2]*vector_2[1];
    c_P[1] = vector_1[2]*vector_2[0] - vector_1[0]*vector_2[2];
    c_P[2] = vector_1[0]*vector_2[1] - vector_1[1]*vector_2[0];
    
    return c_P;
}

double acceleration(double r, double grav_scalar, double centri) { //acceleration function, from newton's first
    accel_new =  -grav_scalar * r;
    return accel_new;
}

double position(double r, double v, double vdot) { //posiiton function 
    r_new = r + v*del_t + 0.5*vdot* pow(del_t, 2.0) ;
    return r_new;
}

double velocity(double v, double vdot, double vdot_next) { //velocity function
    vel_new =  v + 0.5*(vdot + vdot_next)* del_t;
    return vel_new;
}

vector <double> LF_solver(vector < vector <double> > pos, vector < vector <double> > vel, vector < vector <double> > accel,\
                          vector <double> x_positions, vector <double> y_positions, vector <double> z_positions) { //leap-frog integration method
    
    for (unsigned int i = 0; i < 1000; i++) {
        
        //new x coordinates
        x_new = position(pos[i][0], vel[i][0], accel[i][0]);
        x_positions.push_back(x_new);
        
        //new y coordinates
        y_new = position(pos[i][1], vel[i][1], accel[i][1]);
        y_positions.push_back(y_new);
        
        
        //new z coordinates
        z_new = position(pos[i][2], vel[i][2], accel[i][2]);
        z_positions.push_back(z_new);
        
        grav = grav_scalar(x_new, y_new, z_new);
        //omega = cross_product(r0, v0) / pow( pow(x0,2.0) + pow(y0, 2.0) + pow(z0, 2.0), 1.0/2.0)
        
        vxdot_next = acceleration(x_new, grav, centrifugal[0]);
        xdot_new = velocity(vel[i][0], accel[i][0], vxdot_next);
        
        vydot_next = acceleration(y_new, grav, centrifugal[1]);
        ydot_new = velocity(vel[i][1], accel[i][1], vydot_next);
        
        vzdot_next = acceleration(z_new, grav, centrifugal[2]);
        zdot_new = velocity(vel[i][2], accel[i][2], vzdot_next);
        
        
        
        //new positions vector
        position_new.push_back(x_new);
        position_new.push_back(y_new);
        position_new.push_back(z_new);
                               
        pos.push_back(position_new);
        position_new.clear();
        
        //new velocity vector
        velocity_new.push_back(xdot_new);
        velocity_new.push_back(ydot_new);
        velocity_new.push_back(zdot_new);
        
        vel.push_back(velocity_new);
        velocity_new.clear();
            
        //new acceleration vector
        acceleration_new.push_back(vxdot_next);
        acceleration_new.push_back(vydot_next);
        acceleration_new.push_back(vzdot_next);
            
        accel.push_back(acceleration_new);
        acceleration_new.clear();
        
        
    }
    
     POSITIONS = pos;
     VELS = vel;
     ACCELS = accel;
    
     X = x_positions;
     Y = y_positions;
     Z = z_positions;
    
    
     return X;
     return Y;
     return Z;
     
     
    
}

void file_creator(vector<double> x, vector<double> y, vector<double> z) {
	ofstream myfile("test_data.txt");
	if (myfile.is_open()) {
		for (int i = 0; i < x.size(); i++) {
			char string[15];
            
			myfile << x[i] << ",";
			myfile << y[i] << ",";
            myfile << z[i] << "\n";
		}
		myfile.close();
	}
    
	else cout << "unable to open file";

}
               
               
//main code

int main() {
   
    
    //semi-major axis in meters
    a = semimajor(Period_days);
    
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
    double ydot0 = -2.*PI;
    double zdot0 = 0.0; 
    
    //initial velocity vector
    v0.push_back(xdot0);
    v0.push_back(ydot0);
    v0.push_back(zdot0);
    
    gravity_0 = grav_scalar(x0, y0, z0);
    
    //centrifugal force 'pre' term to be multiplied with position in dimensionless units
    //omega_0 = cross_product(r0, v0) / pow( pow(x0,2.0) + pow(y0, 2.0) + pow(z0, 2.0), 1.0/2.0);
    
    centrifugal.push_back(-4.0*pow(PI, 2.0));
    centrifugal.push_back(-4.0*pow(PI, 2.0));
    centrifugal.push_back(0.0);
    
    a0.push_back(acceleration(x0, gravity_0, centrifugal[0]));
    a0.push_back(acceleration(y0, gravity_0, centrifugal[1]));
    a0.push_back(acceleration(z0, gravity_0, centrifugal[2]));
    
    
    //Start 
    positions.push_back(r0);
    velocities.push_back(v0);
    accelerations.push_back(a0);
    
    x_positions.push_back(x0);
    y_positions.push_back(y0);
    z_positions.push_back(z0);
   
    LF_solver(positions, velocities, accelerations, x_positions, y_positions, z_positions);
    
    
    X_positions = X;
    Y_positions = Y;
    Z_positions = Z;
    
    
    
    file_creator(X_positions, Y_positions, Z_positions);
    
    

   
}
