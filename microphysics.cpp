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
  //cout << beta << endl;
  //cout << (4.0/3.0) * k * rho_d * 1.0e-4 << endl; 
	return beta;

}

double clausius_clap(double s, double x, double y, double z, double tau, double Td){
  
  //double T;
  //T_ = brent(s,x,y,z,tau);
  //cout << (exp((-A/Td) + Bp) * pow(amu/(2.0*PI*kb*Td), 0.5) ) << endl;
  //cout << Td << endl;
  double c1 = -alpha*exp((-A/Td) + Bp)/rho_d;
  double c2 = pow((mu*amu)/(2.0*PI*kb*Td), 0.5);
  //cout << Td << endl;
  //cout << "ds1 " << c1 << endl;
  //cout << "ds2 " << c2 << endl;
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

double f_Tdust( double s, double x, double y, double z, double tau, double Tdust){
  double fTdust, dl, sa;
  
  dl = scalar((x-star_pos[0]), y, z)*a*pow(10.0,2);
  //sa = 2.0*(1-pow(1-pow((Rstar*Rsun_cgs)/dl, 2.0), 0.5));
  fTdust = opac.stellar_abs(s) * exp(-tau)*pow(Temp,4) * pow(Rstar*Rsun_cgs, 2)/(4.0*pow(dl,2)) - 
            opac.particle_abs(s, Tdust) * pow(Tdust,4) + opac.particle_abs(s, 100.0) * pow(100.0,4); 

  //cout << opac.particle_abs(Tdust,s) << endl;
  //cout << opac.stellar_abs(s) << endl;
  //fTdust = fTdust/(a*pow(10.0,2));
  return fTdust;
}

double brent(double size, double x, double y, double z, double tau){
  const int itmax=100;
  const double eps = std::numeric_limits<double>::epsilon();
  bool cond1, cond2, cond3, cond4, cond5;
  int mflag=0;
  double Ac, Bc, Cc, Dc;
  double p, q, r, s, tol1, xm;
  double fa, fb, fc, fs;
  double tol = 1.0e-5;
  Ac = 1.0;
  Bc = 5000.0;
  fa = f_Tdust(size,x,y,z,tau, Ac);
  fb = f_Tdust(size,x,y,z,tau, Bc);
  

  if (fa*fb >= 0.0) {
    cout << "root must be bracketed in zbrent" << endl;
    cout << "size " << size << endl;
    cout << "tau " << tau << endl;
    cout << "x " << x << " y " << y << " z " << z << endl;
    cout << "fa " << fa << endl;
    cout << "fb " << fb << endl;
    abort();
  }
  if (abs(fa)<abs(fb)) {
    swap(Ac,Bc);
    swap(fa,fb);
  }
  Cc=Ac;
  fc = fa;
  mflag = 1;
  int iter = 0;
  
  while (abs(Bc-Ac) > tol) {
        if (iter > 100){
          cout << "more than 100 iterations in brent " << endl;
          abort();
        }
        //cout << abs(Bc-Ac) << endl;
        if ((fa != fc) && (fb != fc)) {
          s = ((Ac*fb*fc) / ((fa-fb)*(fa-fc))) +
              ((Bc*fa*fc) / ((fb-fa)*(fb-fc))) +
              ((Cc*fa*fb) / ((fc-fa)*(fc-fb)));
        } else {
          s = Bc - fb*(Bc-Ac)/(fb-fa);
        }
        cond1 = (s < (3.0*Ac + Bc)/4.) || (s > Bc) ;
        cond2 = (mflag==1) && (abs(s-Bc) >= abs(Bc-Cc)/2.);
        cond3 = (mflag==0) && (abs(s-Bc) >= abs(Cc-Dc)/2.);
        cond4 = (mflag==1) && (abs(Bc-Cc) < tol);
        cond5 = (mflag==0) && (abs(Cc-Dc) < tol);

        if (cond1 || cond2 || cond3 || cond4 || cond5 ) {
          s = (Ac + Bc)/2.;
          mflag = 1;
        } else {
          mflag = 0;
        }
        fs = f_Tdust(size,x,y,z,tau,s);
        Dc = Cc;
        Cc = Bc;

        if (fa*fs < 0.0) {
          Bc = s;
        } else {
          Ac = s;
        }

        if (abs(fa)<abs(fb)) {
          swap(Ac, Bc);
          swap(fa, fb);
        }

    iter +=1;
    }
    if (fb == 0.0) {
      return Bc;
    } else {
      return s;
    }
  
  }



double dust_mass(double s){
    double md;
    md = rho_d * (4.0/3.0) * PI * pow(s, 3.0);
    return md;
}

