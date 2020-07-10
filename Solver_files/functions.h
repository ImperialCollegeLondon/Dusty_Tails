//Functions declared here
#include <iostream>
#include <vector>

using namespace std;

///k values declared

extern vector <double> k1, k2, k3, k4, k5, k6, k7;
extern vector <double> k1d, k2d, k3d, k4d, k5d, k6d, k7d;
extern double ks1, ks2, ks3, ks4, ks5, ks6, ks7;



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

//microphysics file
double omega(double mplanet, double mstar);
double beta_fn(double k, double L_star, double M_star);
double opacity(double s, double x, double y, double z);
double qfactor(double s, double x, double y, double z);
double clausius_clap(double s, double x, double y, double z);
double luminosity(double R_star);
double radial_vel(vector <double> vel, vector <double> s_vector);
double temp_dust( double lum, double s,  double x, double y, double z);
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

double acceleration( int i, double pos_star, double pos_planet, vector <double> V);
double sublimation(double s, double x, double y, double z);

vector <double> new_variables(double h, vector <double> V, bool order5);

vector <double> RK_solver(vector <double> V_0, double t_0, double del_t, \
  double h_p);
