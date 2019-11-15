#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "RK_variables.h"

using namespace std;

//some declaration of variables

vector <double> V0; //initial variables vector

vector < vector <double> > Vs, final_variables;

vector <double> x_positions, y_positions, z_positions;

vector <double> x, y, z, V_new;

vector <double> deltas, delta_values;
vector <double> time_plot, semis, t, timing, as;




double semimajor(double period) { //function to evaluate semi-major axis value, period in days, a in AU
    
   semi = pow(7.496e-6 * Mstar_sun * pow(Period_days, 2.0), (1.0/3.0)); //in AU
   return semi*AU_to_m; //output in meters
}

double acceleration( double pos, double x, double y, double z){
     
  vel_dot = (-G_dim * pos) / pow( ( pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) ) , (3.0/2.0) );
  return vel_dot;
}


void k_values(double h, vector <double> V, bool order5){
    
    //k1 values
    k1_xdot = h* acceleration( V[0], V[0], V[1], V[2]);
    k1_x = h* V[3];
    
    k1_ydot = h* acceleration( V[1], V[0], V[1], V[2]);
    k1_y = h* V[4];
    
    k1_zdot = h* acceleration( V[2], V[0], V[1], V[2]);
    k1_z = h* V[5];
    
    //k2 values
    k2_xdot = h*acceleration( V[0] + a21*k1_x, V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z);
    
    k2_x = h* (V[3] + a21*k1_xdot);
    
    k2_ydot = h*acceleration( V[1] + a21*k1_y, V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z);
    
    k2_y = h* (V[4] + a21*k1_ydot);
    
    k2_zdot = h*acceleration( V[2] + a21*k1_z, V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z);
    
    k2_z = h* (V[5] + a21*k1_zdot);
    
    
    //k3 values 
    k3_xdot = h* acceleration( V[0] + a31*k1_x + a32*k2_x, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z);
    
    k3_x = h*(V[3] + a31*k1_xdot + a32*k2_xdot);
    
    k3_ydot = h* acceleration( V[1] + a31*k1_y + a32*k2_y, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z);
    
    k3_y = h*(V[4] + a31*k1_ydot + a32*k2_ydot);
    
    k3_zdot = h* acceleration( V[2] + a31*k1_z + a32*k2_z, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z);
    
    k3_z = h*(V[5] + a31*k1_zdot + a32*k2_zdot);
    
    
    //k4 values
    k4_xdot = h* acceleration( V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z); 
    
    k4_x = h* (V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot);
    
    k4_ydot = h* acceleration( V[1] + a41*k1_y + a42*k2_y + a43*k3_y, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z);
    
    k4_y = h* (V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot);
    
    k4_zdot = h* acceleration( V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z);
    
    k4_z = h* (V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);
    
    //k5 values
    
    k5_xdot = h* acceleration( V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z);
    
    k5_x = h* (V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot);
    
    k5_ydot = h* acceleration( V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z);
    
    k5_y = h* (V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot);
    
    k5_zdot = h* acceleration( V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z);
    
    k5_z = h* (V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);
    
    
    //k6 values
    
    k6_xdot = h* acceleration( V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z);
    
    k6_x = h* (V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot);
    
    k6_ydot = h* acceleration( V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z);
    
    k6_y = h* (V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot);
    
    k6_zdot = h* acceleration( V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z);
    
    k6_z = h* (V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);
    
    
    if (order5 == true) {
        
        k7_xdot = h* acceleration( V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z);
        
        k7_x = h* (V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot);
        
        k7_ydot = h* acceleration( V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z);
        
        k7_y = h* (V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot);
        
        k7_zdot = h* acceleration( V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z);
        
        
        k7_z = h* (V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);
        
    }
    
    
}


double new_variables(double h, vector <double> V, bool order5){
    
    if (order5 == false) {
        
     k_values(h, V, false);
        
     x_new = V[0] + b1*k1_x + b3* k3_x + b4* k4_x + b5* k5_x + b6* k6_x;
     y_new = V[1] + b1*k1_y + b3* k3_y + b4* k4_y + b5* k5_y + b6* k6_y;
     z_new = V[2] + b1*k1_z + b3* k3_z + b4* k4_z + b5* k5_z + b6* k6_z; 
        
     xdot_new = V[3] + b1*k1_xdot  + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;
     ydot_new = V[4] + b1*k1_ydot  + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
     zdot_new = V[5] + b1*k1_zdot  + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;
        
    } else {
        
        k_values(h, V, true);
        
        x_new = V[0] + bs1*k1_x + bs3*k3_x + bs4*k4_x + bs5*k5_x + bs6*k6_x + bs7*k7_x;
        y_new = V[1] + bs1*k1_y + bs3*k3_y + bs4*k4_y + bs5*k5_y + bs6*k6_y + bs7*k7_y;
        z_new = V[2] + bs1*k1_z + bs3*k3_z + bs4*k4_z + bs5*k5_z + bs6*k6_z + bs7*k7_z;
        
        xdot_new = V[3] + bs1*k1_xdot  + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7*k7_xdot;
        ydot_new = V[4] + bs1*k1_ydot  + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7*k7_ydot;
        zdot_new = V[5] + bs1*k1_zdot  + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7*k7_zdot;
        
        
    }
    
    return x_new, y_new, z_new, xdot_new, ydot_new, zdot_new;
}



double delta( double value1, double value2){
        
    del = fabs(value1 - value2) + 1e-10;
   
    return del/tol;
    
      
}

double h_optimal(vector <double> deltas_list, double h){
    
    max_delta = *max_element(deltas_list.begin(), deltas_list.end()); //includes tolerance
    
    h_opt = S * h * pow(1/max_delta, 1.0/5.0);
    
    if (h_opt < h){
        
        return h_opt;
        
    } else return h;

}

vector <double> h_check(double h, vector <double> V){
    
    deltas.clear();
    
    //4th order 
    new_variables(h, V, false);
    
    x4 = x_new;
    y4 = y_new;
    z4 = z_new;
    xdot4 = xdot_new;
    ydot4 = ydot_new;
    zdot4 = zdot_new;

    //5th order 
    new_variables(h, V, true);
    
    x5 = x_new;
    y5 = y_new;
    z5 = z_new;
    xdot5 = xdot_new;
    ydot5 = ydot_new;
    zdot5 = zdot_new;

    
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


vector <double> next_step(double h, vector <double> V) {
    
    //ensure vector are clear
    V_new.clear();
    
    new_variables(h, V, false);
    
    V_new.push_back(x_new);
    V_new.push_back(y_new);
    V_new.push_back(z_new);
    
    V_new.push_back(xdot_new);
    V_new.push_back(ydot_new);
    V_new.push_back(zdot_new);
       

         
    return V_new;
}


void RK_solver(double h0, vector < vector <double> > Vs, \
     vector <double> x_positions, vector <double> y_positions, vector <double> z_positions, vector <double> t){
    
    
    //obtain delta values for the 6 variables
    delta_values = h_check(h0, Vs[0]);
        
        
    //calculate new h
    h_new = h_optimal(delta_values, h0);
        
    next_step(h_new, Vs[0]);
        
    x_positions.push_back(V_new[0]);
    y_positions.push_back(V_new[1]);
    z_positions.push_back(V_new[2]);
        
    Vs.push_back(V_new);
    
    t.push_back(h_new);
        
    
    for (unsigned int i = 1; i < 1e+6; i++) {
        
        //obtain delta values for the 6 variables
        delta_values = h_check(h_new, Vs[i]);
        
        //calculate new h
        h_new = h_optimal(delta_values, h_new);
        t.push_back(t[i] + h_new);
        //cout << h_new << endl;
        
        next_step(h_new, Vs[i]);
        
        x_positions.push_back(V_new[0]);
        y_positions.push_back(V_new[1]);
        z_positions.push_back(V_new[2]);
         
        Vs.push_back(V_new);
       
     } 
    
    
     final_variables = Vs;
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
	ofstream myfile("semi_data.txt");
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

int main() {
    
    
    a = semimajor(Period_days);
    G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant
    
    h0 = 0.001;
    
    //Define initial position in dimensionless units
    
    double x0 = 1.0;
    double y0 = 0.0;
    double z0 = 0.0;
    
    //Define initial velocity in dimensionless units
    
    double xdot0 = 0.0;
    double ydot0 = 2.0*PI;
    double zdot0 = 0.0; 
    
    //initial variables vector
    
    V0.push_back(x0);
    V0.push_back(y0);
    V0.push_back(z0);
    V0.push_back(xdot0);
    V0.push_back(ydot0);
    V0.push_back(zdot0);
    
    Vs.push_back(V0);
    
    x_positions.push_back(x0);
    y_positions.push_back(y0);
    z_positions.push_back(z0);
    
    double t0 = 0.0;
    
    t.push_back(t0);

    RK_solver(h0, Vs, x_positions, y_positions, z_positions, t);
    
    as = semimajor_plot(x, y, z);
    
    file_plot(timing, as);
    
    file_creator(x, y, z);
       
    
}