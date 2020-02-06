#include "Non_uniform_grid_3D_pointers.cpp"


int main() {
  semimajor(Period_days);
  //cout << "Semi-major axis = " << a << endl;
  //cout << "start" << endl;
  //setting constants
  NR=50;
  NP=50;
  NT=50;
  //nr=(NR-1.)/2.;

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
  Ta = new double[51];
  delete[] Tb;
  Tb = new double[51];
  delete[] dTa;
  dTa = new double[51];
  delete[] dTb;
  dTb = new double[51];
  delete[] Pa;
  Pa = new double[51];
  delete[] Pb;
  Pb = new double[51];
  delete[] dPa;
  dPa = new double[51];
  delete[] dPb;
  dPb = new double[51];
  delete[] Ra;
  Ra = new double[51];
  delete[] Rb;
  Rb = new double[51];
  delete[] dRa;
  dRa = new double[51];
  delete[] dRb;
  dRb = new double[51];

  delete[] d;
  d = new double**[51];
  for(int i =0; i<51; i++){
    d[i] = new double*[51];
    for(int j =0; j<51; j++){
      d[i][j] = new double[51];
      for(int k = 0; k<51; k++){
        d[i][j][k] = 0.;
      }
    }
  }

  delete[] kappa;
  kappa = new double**[51];
  for(int i =0; i<51; i++){
    kappa[i] = new double*[51];
    for(int j =0; j<51; j++){
      kappa[i][j] = new double[51];
      for(int k = 0; k<51; k++){
        kappa[i][j][k] = 0.;
      }
    }
  }

  delete[] t;
  t = new double**[51];
  for(int i =0; i<51; i++){
    t[i] = new double*[51];
    for(int j =0; j<51; j++){
      t[i][j] = new double[51];
      for(int k = 0; k<51; k++){
        t[i][j][k] = 0.;
      }
    }
  }

  delete[] T_g;
  T_g = new double[51];
  delete[] P_g;
  P_g = new double[51];
  delete[] R_g;
  R_g = new double[51];
  delete[] DT;
  DT = new double[51];
  delete[] DP;
  DP = new double[51];
  delete[] DR;
  DR = new double[51];
  delete[] T_x;
  T_x = new double[51];
  delete[] P_x;
  P_x = new double[51];
  delete[] R_x;
  R_x = new double[51];
  delete[] T_inv;
  T_inv = new double[51];
  delete[] P_inv;
  P_inv = new double[51];
  delete[] R_inv;
  R_inv = new double[51];

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
  file_creator_t(t);
  // //file_creator_gauss(gauss,Rb,NR);
  // file_creator_DR(DR,Ra);
  //file_creator_xRa(x,Ra);
  file_creator_d(d,Ta,Pa,Ra,NT,NP,NR);
  //cout << "BAA5" << endl;

  return 0;
}
