//variables for RK solver declared here
#include <vector>

double a, G_dim, h0, t0, semi, h_new, vel_dot; // semi-major axis, dimensionless G and initial time step

//k values

double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x;
double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y;
double k1_z, k2_z, k3_z, k4_z, k5_z, k6_z, k7_z;

double k1_xdot, k2_xdot, k3_xdot, k4_xdot, k5_xdot, k6_xdot, k7_xdot;
double k1_ydot, k2_ydot, k3_ydot, k4_ydot, k5_ydot, k6_ydot, k7_ydot;
double k1_zdot, k2_zdot, k3_zdot, k4_zdot, k5_zdot, k6_zdot, k7_zdot;

double xdot_new, ydot_new, zdot_new;
double x_new, y_new, z_new;

double xdot4, ydot4, zdot4, x4, y4, z4;
double xdot5, ydot5, zdot5, x5, y5, z5;

double xdot_err, ydot_err, zdot_err;
double x_err, y_err, z_err;

double del, h_opt, max_delta, scalar;

