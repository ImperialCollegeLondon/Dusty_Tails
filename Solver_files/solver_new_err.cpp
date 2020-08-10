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
  drag = pr_drag(V)[i];

  vel_dot = grav_star - centri - 2.0*coriol  + radiation - drag + grav_planet;

  return vel_dot;
}

double sublimation(double s, double x, double y, double z){
  double ds1, ds2, ds;

  ds1 = ((-alpha*clausius_clap(s, x, y, z))/rho_d);

  ds2 = pow((mu*amu)/(2.0*PI*kb*temp_dust(luminosity(Rstar), s, x, y, z)), 0.5);
  ds = ds1 * ds2;

  return ds;

}

vector <double> new_variables(double h, vector <double> V, bool order5){
    //function to obtain next step values
    vector <double> new_pos(3), new_vel(3), new_vars(7);
    vector <double> V0 = {V[0], V[1], V[2]};
    vector <double> V0dot = {V[3], V[4], V[5]};
    double s_new;

    new_vars.clear();

    if (order5 == false) {
        k_values(h, V, false, k1, k2, k3, k4, k5, k6, k7, k1d, k2d, k3d, k4d, k5d, k6d, k7d);

        for (unsigned int i = 0; i < 3; i++) {

         new_pos[i] = V0[i] + b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i];
         new_vel[i] = V0dot[i] + b1*k1d[i] + b3*k3d[i] + b4*k4d[i] + b5*k5d[i] + b6*k6d[i];

        }

        s_new = V[6] + b1*ks1 + b3*ks3 + b4*ks4 + b5*ks5 + b6*ks6;

    } else {
        k_values(h, V, true, k1, k2 , k3, k4, k5 , k6, k7, k1d, k2d, k3d, k4d, k5d, k6d, k7d);

        for (unsigned int i = 0; i < 3; i++) {
            new_pos[i] = V0[i] + bs1*k1[i] + bs3*k3[i] + bs4*k4[i] + bs5*k5[i] + bs6*k6[i] + bs7*k7[i];
            new_vel[i] = V0dot[i] + bs1*k1d[i] + bs3*k3d[i] + bs4*k4d[i] + bs5*k5d[i] + bs6*k6d[i] + bs7*k7d[i];
        }
        s_new = V[6] + bs1*ks1 + bs3*ks3 + bs4*ks4 + bs5*ks5 + bs6*ks6 + bs7*ks7;
    }
    new_vars = {new_pos[0], new_pos[1], new_pos[2], new_vel[0], new_vel[1], new_vel[2], s_new};
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

    double t = t_0;
    double maximum_err;
    vector <double> step_sizes;
    double old_h, new_h;
    vector <double> new_vector;
    new_vector.clear();
    step_sizes.clear();

    //obtain delta values for the 6 variables
    maximum_err = error_max(h_p, V_0);

    step_sizes = new_step_size(maximum_err, h_p, 0, V_0);
    //first value is old step size, second value is step size to be used in the next iteration
    old_h = step_sizes[0];
    new_h = step_sizes[1];

    if (new_h == -1.0){
      
      t = del_t;
      return {V_0[0], V_0[1], V_0[2], V_0[3], V_0[4], V_0[5], 0.01e-4, new_h};
    }

    new_vector = next_step(old_h, V_0);
    t = t + old_h;

    if (new_vector[6] < 0.1e-4){
      
      t = del_t;
      return {new_vector[0], new_vector[1], new_vector[2], new_vector[3], new_vector[4], new_vector[5], 0.01e-4, new_h};
    }


    if (t > del_t){
      new_vector = next_step(del_t - (t-old_h), new_vector);
      new_vector.push_back(new_h);
      return new_vector;
    }

    while (t < del_t) {

        //obtain delta values for the 6 variables

        maximum_err = error_max(new_h, new_vector);
        step_sizes = new_step_size(maximum_err, new_h, 0, new_vector);
        old_h = step_sizes[0];
        new_h = step_sizes[1];
        if (new_h == -1.0){
         
          t = del_t;
          return {new_vector[0], new_vector[1], new_vector[2], new_vector[3], new_vector[4], new_vector[5], 0.01e-4, new_h};
        }
        t = t + old_h;
        if (t < del_t) {
             new_vector = next_step(old_h, new_vector);

             if (new_vector[6] < 0.1e-4){
               
               t = del_t;
               return {new_vector[0], new_vector[1], new_vector[2], new_vector[3], new_vector[4], new_vector[5], 0.01e-4, new_h};
             }
           } else {
             new_vector = next_step(del_t - (t-old_h), new_vector);

             if (new_vector[6] < 0.1e-4){
               
               return {new_vector[0], new_vector[1], new_vector[2], new_vector[3], new_vector[4], new_vector[5], 0.01e-4, new_h};
             } else {
             
             new_vector.push_back(new_h);
             return new_vector;
           }
        }
    }
}