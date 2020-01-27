#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"


using namespace std;


double acceleration( double pos_star, double pos_planet, double x, double y, double z, \
                     double centri, double coriol, double radiation, double drag){

  double vel_dot;

  double grav_star, grav_planet;

  grav_star = (-G_dim * pos_star ) / pow( scalar(x - star_x, y ,z), 3.0 );

  grav_planet = (-G_dim * m_planet* pos_planet ) / pow( scalar(x - planet_x + r_planet_dim,y ,z), 3.0 );

  vel_dot = grav_star - centri - 2.0*coriol - drag + radiation + grav_planet;

  return vel_dot;
}

vector <double> new_variables(double h, vector <double> V, bool order5){
    //function to obtain next step values

    vector <double> new_vars;
    double x_new, y_new, z_new;
    double xdot_new, ydot_new, zdot_new;

    new_vars.clear();

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

    new_vars.push_back(x_new);
    new_vars.push_back(y_new);
    new_vars.push_back(z_new);
    new_vars.push_back(xdot_new);
    new_vars.push_back(ydot_new);
    new_vars.push_back(zdot_new);

    return new_vars;
}


vector <double> next_step(double h, vector <double> V){
    //ensure vector is  clear
    vector <double> V_new;
    V_new.clear();

    V_new = new_variables(h, V, false);

    return V_new;
}


vector <double> RK_solver(vector <double> V_0, double t_0, \
               double del_t, double h_p){

    double h_old, h_new;
    double t = t_0;


    vector <double> new_vector, delta_values;

    //obtain delta values for the 6 variables
    delta_values = h_check(h_p, V_0);

    h_old = h_p;

    //calculate new h
    h_new = h_optimal(delta_values, h_old);

    if (h_new < 0.0){
       new_vector = next_step(h_old, V_0);
       t = t + h_old;
       h_new = 2.0*h_old;

     } else {
       new_vector = next_step(h_new, V_0);
       t = t + h_new;
       h_new = h_new;
     }


    while (t < del_t) {

        //cout << t << endl;
        //obtain delta values for the 6 variables
        delta_values = h_check(h_new, new_vector);
        h_old = h_new;

        //calculate new h
        h_new = h_optimal(delta_values, h_new);

        if (h_new < 0.0){
           t = t + h_old;
           if (t < del_t) {
             new_vector = next_step(h_old, new_vector);
             h_new = 2.*h_old;
           } else {
             new_vector = next_step(del_t - (t-h_old), new_vector);
             return new_vector;
           }

         } else {
           t = t + h_new;
           if (t < del_t) {
             new_vector = next_step(h_new, new_vector);
             h_new = h_new;
           } else {
              new_vector = next_step(del_t - (t-h_new), new_vector);
              return new_vector;

           }

     }

  }
}
