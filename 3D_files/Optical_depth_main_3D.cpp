#include "file_creators.cpp"


int main() {
  semimajor(Period_days);
  //cout << "Semi-major axis = " << a << endl;

  //cout << "start" << endl;
  //setting constants
  NR=200;
  NP=200;
  NT=200;
  //nr=(NR-1.)/2.;

  rmin = 0.;
  rmax = 2.*a;
  Rmin = rmin/a;
  Rmax = rmax/a;

  Pmin = 0.;
  Pmax = 2*M_PI;

  Tmin = 0.;
  Tmax = M_PI;

  R_B=0.9; //A and B are in terms of a
  P_B=0.9;
  T_B=0.9;

  R_C=6.; //sd of gaussian grid = 2/C
  P_C=6.;
  T_C=6.;

  R_mu=1.;
  R_sd=Rmax/R_C;
  P_mu=1.;
  P_sd=Pmax/P_C;
  T_mu=1.;
  T_sd=Tmax/T_C;

  //for density
  mean = 1.; //in dimensionless units - this should be at a.
  stde = 0.2*mean; //arbitrary

  //clearing vectors
  //cout << "delete start" << endl;
  delete[] Ta;
  Ta = new double[NT+1];
  delete[] Tb;
  Tb = new double[NT+1];
  delete[] dTa;
  dTa = new double[NT+1];
  delete[] dTb;
  dTb = new double[NT+1];
  delete[] Pa;
  Pa = new double[NP+1];
  delete[] Pb;
  Pb = new double[NP+1];
  delete[] dPa;
  dPa = new double[NP+1];
  delete[] dPb;
  dPb = new double[NP+1];
  delete[] Ra;
  Ra = new double[NR+1];
  delete[] Rb;
  Rb = new double[NR+1];
  delete[] dRa;
  dRa = new double[NR+1];
  delete[] dRb;
  dRb = new double[NR+1];

  delete[] d;
  d = new double**[NT+1];
  for(int i = 0; i<(NT+1); i++){
    d[i] = new double*[NP+1];
    for(int j = 0; j<(NP+1); j++){
      d[i][j] = new double[NR+1];
      for(int k = 0; k<(NR+1); k++){
        d[i][j][k] = 0.;
      }
    }
  }

  delete[] kappa;
  kappa = new double**[NT+1];
  for(int i =0; i<(NT+1); i++){
    kappa[i] = new double*[NP+1];
    for(int j =0; j<(NP+1); j++){
      kappa[i][j] = new double[NR+1];
      for(int k = 0; k<(NR+1); k++){
        kappa[i][j][k] = 0.;
      }
    }
  }

  delete[] t;
  t = new double**[NT+1];
  for(int i =0; i<(NT+1); i++){
    t[i] = new double*[NP+1];
    for(int j =0; j<(NP+1); j++){
      t[i][j] = new double[NR+1];
      for(int k = 0; k<(NR+1); k++){
        t[i][j][k] = 0.;
      }
    }
  }

  delete[] T_g;
  T_g = new double[NT+1];
  delete[] P_g;
  P_g = new double[NP+1];
  delete[] R_g;
  R_g = new double[NR+1];
  delete[] DT;
  DT = new double[NT+1];
  delete[] DP;
  DP = new double[NP+1];
  delete[] DR;
  DR = new double[NR+1];
  delete[] T_x;
  T_x = new double[NT+1];
  delete[] P_x;
  P_x = new double[NP+1];
  delete[] R_x;
  R_x = new double[NR+1];
  delete[] T_inv;
  T_inv = new double[NT+1];
  delete[] P_inv;
  P_inv = new double[NP+1];
  delete[] R_inv;
  R_inv = new double[NR+1];

  //cout << "delete end" << endl;

  //adding initial conditions
  //cout << "BEE" << endl;
  Ta[0] = 0.;
  dTb[0] = 1.;
  Pa[0] = 0.; //segmentation fault here ??
  dPb[0] = 1.;
  Ra[0] = 0.;
  dRb[0] = 1.;
  //t[0] = 0.;

  //running code
  //cout << "starting" << endl;
  //cout << sd << endl;
  make_gauss(NT, NP, NR);
  //cout << "BAAA" << endl;
  //cout << "done one" << endl;
  find_inv(NT, NP, NR,T_B, P_B, R_B, T_g, P_g, R_g);
  //cout << "BAA1" << endl;
  find_DR(NT, NP, NR, T_A=(2./T_suminv), P_A=(2./P_suminv), R_A=(2./R_suminv), T_inv, P_inv, R_inv);
  //cout << "BAA2" << endl;
  // file_creator_gaussian(g,x);
  // file_creator_DR(DR,x);
  build_grid(NT, NP, NR, DT, DP, DR);
  //cout << "BAA3" << endl;
  density_fill(NT, NP, NR,mean,stde);
  //cout << "BAA4" << endl;
  opacity_fill(NT, NP, NR);
  calculate_optical_depth(NT, NP, NR, kappa, d);

  //file creators
  file_creator_t(NT,NP,NR,t);
  // //file_creator_gauss(gauss,Rb,NR);
  // file_creator_DR(DR,Ra);
  //file_creator_xRa(x,Ra);
  file_creator_d(d,Ta,Pa,Ra,NT,NP,NR);
  //cout << "BAA5" << endl;

  return 0;
}
