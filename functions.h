//Functions declared here
#include <iostream>
#include <vector>
#include "spline.h"

using namespace std;

///k values declared



extern vector <double> k1, k2, k3, k4, k5, k6, k7;
extern vector <double> k1d, k2d, k3d, k4d, k5d, k6d, k7d;
extern double ks1, ks2, ks3, ks4, ks5, ks6, ks7;

extern vector < vector < tk:: spline > > s_phi;


extern double r_min, r_max, theta_min, theta_max, phi_min, phi_max;

extern double dr, dtheta, dphi;
//extern double d_r_min, d_r_max, d_t_max, d_t_min, d_p_min, d_p_max, d_dr, d_dtheta, d_dphi;

extern vector < vector < vector <double> > >  tau;
extern vector <double> radii_v, thetas_v, phis_v;


double fRand(double fMin, double fMax);

//forces file
vector <double> centrifugal(vector <double> V);

vector <double> coriolis(vector <double> V);

vector <double> rad_pressure(vector <double> V);

vector <double> pr_drag(vector <double> V);


//maths file

double scalar(double x, double y, double z);

vector <double> cross_product(double m1, double m2, double m3, \
                              double n1, double n2, double n3);

double dot_product(vector <double> n,  vector <double> m);
vector <double> pos_to_spherical(double x, double y, double z);
vector <double> vel_to_spherical(double x, double y, double z);

//microphysics file
double omega(double mplanet, double mstar);
double beta_fn(double k, double tau);
double opacity(double s, double x, double y, double z);
double qfactor(double s, double x, double y, double z);
double clausius_clap(double s, double x, double y, double z);
//double luminosity(double R_star);
double radial_vel(vector <double> vel, vector <double> s_vector);
double temp_dust( double s,  double x, double y, double z);
double dust_mass( double s);

vector <double> drag_vel(vector <double> V);
vector <double> sunit_vector(vector <double> V);

//kvalues file

void k_values(double h, vector <double> V, bool order5, vector <double> &k1, \
    vector <double> &k2, vector <double> &k3, vector <double> &k4, vector <double> &k5, \
    vector <double> &k6, vector <double> &k7, vector <double> &k1d, \
    vector <double> &k2d, vector <double> &k3d, vector <double> &k4d, vector <double> &k5d, \
    vector <double> &k6d, vector <double> &k7d);

//errors

double error( double value1, double value2);
double error_max(double h, vector <double> V);
vector <double> new_step_size(double max_err, double h_old, int fail_status, vector <double> V);




//solver

double acceleration( int i, double pos_star, double pos_planet, vector <double> V, \
                      double centrif, double coriolis, double rad_press, double dragf);
double sublimation(double s, double x, double y, double z);

vector <double> new_variables(double h, vector <double> V, bool order5);

vector <double> RK_solver(vector <double> V_0, double t_0, double del_t, \
  double h_p);

//ray tracer

void build_grids(double *r_a_grid, double *r_b_grid, double *theta_a_grid, \
                double *theta_b_grid, double dr, double dtheta, double dphi,
                double *phi_a_grid, double *phi_b_grid, double r_start, double theta_start, double phi_start);

double gauss(double var, double var_c, double std);



vector <double> grid_scaling(vector <double> s_position);

double r_reverse(double old_r);

double theta_reverse(double old_theta);

double phi_reverse(double old_phi);


double cubicInterpolate ( vector <double> p, vector <double> s1, double x);

double bicubicInterpolate (vector< vector <double> > p,vector <double> s1, vector <double> s2, double x, double y);

double tricubicInterpolate (vector <vector< vector <double> > > p, vector <double> s1, vector <double> s2, vector <double> s3, double x, double y, double z);

double tau_p (vector<double> theta_splines_p, vector<double> radii, vector<double> thetas, vector<double> phis, double radius);
vector < vector < tk:: spline > > splines_phi (vector <vector< vector <double> > > taus, vector <double> radii, vector <double> thetas, vector <double> phis);
vector <vector <double> > phi_spline_result(vector < vector < tk:: spline >> splines , vector<double> radii, \
                            vector <double> thetas, vector<double> phis, double phi);
vector < tk:: spline > splines_theta( vector <vector <double>> phi_splines_p, vector<double> radii, vector<double> thetas);
vector <double> theta_spline_result( vector <tk:: spline> splines, vector<double> radii, vector<double> thetas, vector<double> phis , double theta);

vector <double> vel_grid_scaling(vector <double> s_velocity);