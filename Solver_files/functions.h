//Functions declared here
#include <iostream>
#include <vector>
using namespace std;


///k values declared
extern double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x;
extern double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y;
extern double k1_z, k2_z, k3_z, k4_z, k5_z, k6_z, k7_z;

extern double k1_xdot, k2_xdot, k3_xdot, k4_xdot, k5_xdot, k6_xdot, k7_xdot;
extern double k1_ydot, k2_ydot, k3_ydot, k4_ydot, k5_ydot, k6_ydot, k7_ydot;
extern double k1_zdot, k2_zdot, k3_zdot, k4_zdot, k5_zdot, k6_zdot, k7_zdot;

//forces file
vector <double> centrifugal(double x, double y, double z, \
                            double vx, double vy, double vz);

vector <double> coriolis(double x, double y, double z, \
                         double vx, double vy, double vz);

vector <double> rad_pressure(double x, double y, double z, \
                             double vx, double vy, double vz);

vector <double> pr_drag(double x, double y, double z, \
                        double vx, double vy, double vz);


//maths file

vector <double> to_vector(double x, double y, double z);

double scalar(double x, double y, double z);

vector <double> cross_product(double m1, double m2, double m3, \
                              double n1, double n2, double n3);

double dot_product(vector <double> n,  vector <double> m);

//microphysics file
double omega(double mplanet, double mstar);
double beta_fn(double k, double L_star, double M_star);
double opacity(double Q_fn);
double luminosity(double R_star);
double radial_vel(vector <double> vel, vector <double> s_vector);

vector <double> drag_vel(double x, double y, double z, double vx, double vy, double vz);
vector <double> sunit_vector(double x, double y, double z);

//kvalues file

void k_values(double h, vector <double> V, bool order5);

//errors

double delta( double value1, double value2);
double h_optimal(vector <double> deltas_list, double h);
vector <double> h_check(double h, vector <double> V);

//solver

double acceleration( double pos_star, double pos_planet, double x, double y, double z, \
                     double centri, double coriol, double radiation, double drag);

vector <double> new_variables(double h, vector <double> V, bool order5);

vector <double> next_step(double h, vector <double> V);

vector <double> RK_solver(vector <double> V_0, double t_0, double del_t, \
  double h_p);
