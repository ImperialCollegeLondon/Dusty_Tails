#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
//#include "RK_variables.h"
#include "functions.h"

using namespace std;

double xdot_new, ydot_new, zdot_new;
double x_new, y_new, z_new;



double acceleration( double pos_star, double pos_planet, double x, double y, double z, \
                     double coriol, double centri, double radiation, double drag, \
                     vector<double> star, vector <double> planet){

  //double Mplanet = (300.0 * Mearth) / (Mstar_kg);
  //function to solve final equation
  double vel_dot;
  vel_dot = ((-G_dim * 1.0 * pos_star ) / pow( scalar(x - star[0] ,y - star[1],z - star[2]), 3.0 )) \
					- centri - 2.0*coriol;
  // + ((-G_dim * (Mplanet)* pos_planet ) / pow( scalar(x - planet[0],y - planet[1],z - planet[2]), 3.0 ))
  return vel_dot;
}

double new_variables(double h, vector <double> V, bool order5, \
                      vector <double> star, vector <double> planet){
    //function to obtain next step values
    if (order5 == false) {

     k_values(h, V, false, star, planet);

     x_new = V[0] + b1*k1_x + b3* k3_x + b4* k4_x + b5* k5_x + b6* k6_x;
     y_new = V[1] + b1*k1_y + b3* k3_y + b4* k4_y + b5* k5_y + b6* k6_y;
     z_new = V[2] + b1*k1_z + b3* k3_z + b4* k4_z + b5* k5_z + b6* k6_z;

     xdot_new = V[3] + b1*k1_xdot  + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;
     ydot_new = V[4] + b1*k1_ydot  + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
     zdot_new = V[5] + b1*k1_zdot  + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;

    } else {

        k_values(h, V, true, star, planet);

        x_new = V[0] + bs1*k1_x + bs3*k3_x + bs4*k4_x + bs5*k5_x + bs6*k6_x + bs7*k7_x;
        y_new = V[1] + bs1*k1_y + bs3*k3_y + bs4*k4_y + bs5*k5_y + bs6*k6_y + bs7*k7_y;
        z_new = V[2] + bs1*k1_z + bs3*k3_z + bs4*k4_z + bs5*k5_z + bs6*k6_z + bs7*k7_z;

        xdot_new = V[3] + bs1*k1_xdot  + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7*k7_xdot;
        ydot_new = V[4] + bs1*k1_ydot  + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7*k7_ydot;
        zdot_new = V[5] + bs1*k1_zdot  + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7*k7_zdot;


    }

    return x_new, y_new, z_new, xdot_new, ydot_new, zdot_new;
}


vector <double> next_step(double h, vector <double> V, vector <double> star, \
                          vector <double> planet) {

    //ensure vector is  clear
    V_new.clear();

    new_variables(h, V, false, star, planet);

    V_new.push_back(x_new);
    V_new.push_back(y_new);
    V_new.push_back(z_new);

    V_new.push_back(xdot_new);
    V_new.push_back(ydot_new);
    V_new.push_back(zdot_new);

    return V_new;
}


void RK_solver(double h0, vector <double> V_0, double t_0, vector <double> star, \
               vector <double> planet){

    double h_old;
    double t = t_0;
    double m_planet;
    double r_h;


    //m_planet = (300.0 * Mearth) / (Mstar_kg); //planet mass in terms of star mass
    m_planet = 0.0;
    r_h = pow(m_planet/3.0, 1.0/3.0); //Hill radius
    ofstream file("planet_data7.txt");

    double planetary_pos = 1.0 / (m_planet + 1.0);
    file << t << ",";
    //file << fabs((0.1*r_h - scalar(V_0[0] - planetary_pos, V_0[1] , V_0[2]))) / (0.1*r_h)<< "\n";
    file << V_0[0] << ",";
    file << V_0[1] << ",";
    file << V_0[2] << "\n";

    //obtain delta values for the 6 variables
    delta_values = h_check(h0, V_0, star, planet);

    h_old = h0;


    //calculate new h
    h_new = h_optimal(delta_values, h_old);

    if (h_new < 0.0){
       next_step(h_old, V_0, star, planet);
       t = t + h_old;
       h_new = h_old;
     } else {
       next_step(h_new, V_0, star, planet);
       t = t + h_new;
       h_new = h_new;
     }

    double next_time = 0.1;

    double delta_time = 0.1;

    while (next_time < 5.1){

        //obtain delta values for the 6 variables
        delta_values = h_check(h_new, V_new, star, planet);
        h_old = h_new;
        //calculate new h
        h_new = h_optimal(delta_values, h_new);

        if (h_new < 0.0){
           next_step(h_old, V_new, star, planet);
           t = t + h_old;
           h_new = 2.*h_old;
           //cout << "old h twice " << 2.0*h_old << endl;

         } else {
           next_step(h_new, V_new, star, planet);
           t = t + h_new;
           h_new = h_new;
           //cout << "new h " << h_new << endl;
         }

        if ( t > next_time) {

            file << t << ",";
            //file << fabs(0.1*r_h - (scalar(V_new[0]- planetary_pos, V_new[1] , V_new[2])))  / (0.1*r_h) << "\n";
	          file << V_new[0] << "," ;
	          file << V_new[1] << ",";
	          file << V_new[2] << "\n";

	          next_time = next_time + delta_time;
        }

     }

}
