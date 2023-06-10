#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include "opacities.h"
#include <omp.h>
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

using namespace std;
//initialisation for tables for opacity calculations


double beta_fn(double k, double tau, double s){
  //mass of star in terms of mass of sun, and everything in cgs units
  //function to evaluate beta (radiation accel/ gravity accel)
  double beta =0.0, qfac;
  if (tau < 1.0e-33) {
    tau = 0.0;
  }
  qfac = (4.0/3.0) * k * rho_d * s;
	beta = (3./16.)*(qfac*lum*exp(-tau)) / (rho_d*s*PI*c_cgs*G_cgs*Mstar_sun*Msun_cgs);

	return beta;

}

double clausius_clap(double s, double x, double y, double z, double tau, double Td, bool debug){
  
  double c1 = -alpha*exp((-A/Td) + Bp)/rho_d;
  double c2 = pow((mu*amu)/(2.0*PI*kb*Td), 0.5)*T;
  if (debug) {
    cout << "c1 " << c1 << endl;
    cout << "c2 " << c2 << endl;
  }
  return c1*c2;
  }

double radial_vel(vector <double> vel, vector <double> s_vector){
  //function to obtain radial velocity term
	return dot_product(vel, s_vector);

}

vector <double> drag_vel(vector <double> V){
   //function to obtain velocity relative to radiation field
	 vector <double> omegar, vrad(3);
	 omegar = cross_product(0.0, 0.0, ang_vel, V[0], V[1], V[2]);

     for (unsigned int i=0; i < 3; i++){
         vrad[i] = V[i+3] + omegar[i];
    }
	 return vrad;
}

vector <double> sunit_vector(vector <double> V){
  vector <double> s_unit(3);

  for (unsigned int i = 0; i < 3; i++) {
  //function to obtain unit vector of radiation field
    s_unit[i] = (V[i] - star_pos[i]) / scalar(V[0] - star_pos[0],V[1],V[2]);

   }

	return s_unit;
}

double f_Tdust( double s, double dl, double tau, double Tdust){
  double fTdust, sa;
  
  fTdust = opac.stellar_abs(s) * exp(-tau)*pow(Temp,4.) * pow(Rstar*Rsun_cgs, 2.)/(4.0*pow(dl,2.)) - 
            opac.particle_abs(s, Tdust) * pow(Tdust,4.) + opac.particle_abs(s, 100.0) * pow(100.0,4.); 

  return fTdust;
}

double brent(double size, double x, double y, double z, double tau){
  //Brent's method base on Numerical Recipes, to solve for the temperature of a dust grain
  const int itmax=100;
  const double eps = 1.0e-8;
  bool cond1, cond2, cond3, cond4, cond5;
  int mflag=0, iter;
  double T_d_guess, dl;
  double Ac, Bc, Cc, Dc, Ec;
  double p, q, r, s, tol1, xm;
  double fa, fb, fc, fs;
  double tol = 1.0e-2;
  double min1, min2;
  dl = scalar((x-star_pos[0]), y, z)*a;

  T_d_guess = Temp * pow(exp(-tau), 0.25) * pow(Rstar*Rsun_cgs, 0.5)/pow(2.0*dl, 0.5);
  //cout << T_d_guess << endl;
  Ac = 100.0;
  Bc = 5000.0;
  Cc=Bc;
  fa = f_Tdust(size,dl,tau, Ac);
  fb = f_Tdust(size,dl,tau, Bc);
  

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    cout << "root must be bracketed in zbrent" << endl;
    abort();
  }
  fc = fb;
  for (iter=1;iter<=itmax;iter++){
  if ((fb > 0.0  && fc > 0.0) || (fb < 0.0 && fc < 0.0)){
    Cc=Ac;
    fc=fa;
    Ec = Dc = Bc-Ac;
  }

  if (abs(fc)<abs(fb)) {
    Ac = Bc;
    Bc = Cc;
    Cc = Ac;
    fa = fb;
    fb = fc;
    fc = fa;
  }
  tol1 = 2.0*eps*abs(Bc)+0.5*tol;
  //cout << "tolerance " << tol1 << endl;
  xm = 0.5*(Cc-Bc);
  //cout << "xm " << xm << endl;
  if (abs(xm) <= tol1 || fb == 0.0){
    //cout << "Cc " << Cc << " Bc " << Bc << endl;
    return Bc;
  }
  if (abs(Ec) >= tol1 && abs(fa)>abs(fb)) {
    //attempt inverse quadratic interpolation
    s = fb/fa;
    if (Ac == Cc) {
      p = 2.0*xm*s;
      q = 1.0-s;
    } else{
      q = fa/fc;
      r = fb/fc;
      p = s*(2.0*xm*q*(q-r)-(Bc-Ac)*(r-1.0));
      q = (q-1.0)*(r-1.0)*(s-1.0);
    }
    if (p>0.0) {
      q = -1.0*q;
    }
    p = abs(p);
    min1 = 3.0*xm*q-abs(tol1*q);
    min2 = abs(Ec*q);

    if (2.0*p < (min1 < min2 ? min1 : min2)) {
      //accept interpolation
      Ec = Dc;
      Dc = p/q;
    } else{
      //interpolation failed, use bisection.
      Dc = xm;
      Ec = Dc;
    }
  } else {
    //bounds decreasing too slowly, use bisection.
    Dc = xm;
    Ec = Dc;
  }
  Ac = Bc;
  fa = fb;
  if (abs(Dc)>tol1) {
    Bc += Dc;
  } else{
    Bc += SIGN(tol1, xm);
  }
  fb = f_Tdust(size,dl,tau, Bc);
  }
  cout << "Maximum number of iterations reached." << endl;
  //cout << "size " << size << endl;
  //cout << "x " << x << " y " << y << " z " << z << endl;
  abort();
  }



double dust_mass(double s){
    double md;
    md = rho_d * (4.0/3.0) * PI * pow(s, 3.0);
    return md;
}

