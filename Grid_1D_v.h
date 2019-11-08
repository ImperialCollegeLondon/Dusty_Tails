// Defining Variables

#include <vector>
#include "constants.h"

using namespace std;

//grid boundaries

double rmin=0.;
double Rmin=rmin/Rsun;
double rmax=2.0*Rsun;
double Rmax=rmax/Rsun;

// grid sizes

double NR=100.; // number of grid cells?

// grid starting/ending indicies

int is;
int ie;
int itot;

// grid coordinates (using vectors as arrays)
vector <double> Ra; // a-mesh defines edges of grid cells
vector <double> Rb; // b-mesh defines centers of grid cells
vector <double> dRa;
vector <double> dRb;
double Ra_new, Rb_new, dRa_new, dRb_new;

// density

vector <double> rho; //density
vector <double> d; //dimensionless density

// opacity

vector <double> kappa; //opacity
vector <double> k; //dimensionless opacity

//optical depth

vector <double> t; //optical depth (dimensionless anyway)
double t_new;
