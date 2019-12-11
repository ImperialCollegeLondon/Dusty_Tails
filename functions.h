//Functions declared here
#include <vector>

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

//NOT DONE YET
