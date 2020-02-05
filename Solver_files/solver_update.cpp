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

double star_x = -((m_planet) / (m_planet + 1.0));
double planet_x = 1.0 / (m_planet + 1.0);

vector <double> star_pos = {star_x,0.0, 0.0 };

vector <double> planet_pos = {planet_x, 0.0, 0.0};


double acceleration( int i, double pos_star, double pos_planet, vector <double> V){

  double vel_dot, grav_star, grav_planet, coriol, drag, radiation, centri;

  grav_star = (-G_dim * pos_star ) / pow( scalar(V[0] - star_pos[0], V[1] , V[2]), 3.0);

  grav_planet = (-G_dim * m_planet* pos_planet ) / pow( scalar(V[0] - planet_pos[0]  + r_planet_dim, V[1] , V[2]), 3.0 );

  coriol = coriolis(V)[i];
  centri = centrifugal(V)[i];
  radiation = rad_pressure(V)[i];
  //cout << "radiation " << radiation << "in " << i << endl;
  drag = pr_drag(V)[i];

  //cout << "grav " << grav_star << endl;
  //cout << "centri " << centri << endl;

  vel_dot = grav_star - centri - 2.0*coriol  + radiation - drag + grav_planet;

  return vel_dot;
}

vector <double> new_variables(double h, vector <double> V, bool order5){
    //function to obtain next step values
    vector <double> new_pos(3), new_vel(3), new_vars(6);
    vector <double> V0 = {V[0], V[1], V[2]};
    vector <double> V0dot = {V[3], V[4], V[5]};

    new_vars.clear();

    if (order5 == false) {
        k_values(h, V, false, k1, k2, k3, k4, k5, k6, k7, k1d, k2d, k3d, k4d, k5d, k6d, k7d);

        for (unsigned int i = 0; i < 3; i++) {

         new_pos[i] = V0[i] + b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i];
         new_vel[i] = V0dot[i] + b1*k1d[i] + b3*k3d[i] + b4*k4d[i] + b5*k5d[i] + b6*k6d[i];
        }
    } else {
        k_values(h, V, true, k1, k2 , k3, k4, k5 , k6, k7, k1d, k2d, k3d, k4d, k5d, k6d, k7d);

        for (unsigned int i = 0; i < 3; i++) {
            new_pos[i] = V0[i] + bs1*k1[i] + bs3*k3[i] + bs4*k4[i] + bs5*k5[i] + bs6*k6[i] + bs7*k7[i];
            new_vel[i] = V0dot[i] + bs1*k1d[i] + bs3*k3d[i] + bs4*k4d[i] + bs5*k5d[i] + bs6*k6d[i] + bs7*k7d[i];
        }
    }
    new_vars = {new_pos[0], new_pos[1], new_pos[2], new_vel[0], new_vel[1], new_vel[2]};
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

    //cout << "h_new " << h_new << endl;

    if (h_new < 0.0){
       new_vector = next_step(h_old, V_0);
       t = t + h_old;
       h_new = h_old;

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
             h_new = h_old;
           } else {
             new_vector = next_step(del_t - (t-h_old), new_vector);
             new_vector.push_back(h_old);
             return new_vector;
           }

         } else {
           t = t + h_new;
           if (t < del_t) {
             new_vector = next_step(h_new, new_vector);
             h_new = h_new;
           } else {
              new_vector = next_step(del_t - (t-h_new), new_vector);
              new_vector.push_back(h_new);
              return new_vector;

           }

     }

  }
}
