#include "file_creators.cpp"


int main() {
  semimajor(Period_days);
  //cout << "Semi-major axis = " << a << endl;

  for(noparticles=1000.;noparticles<=10000.;(noparticles=noparticles+1000.)){

    //cout << "start" << endl;
    //setting constants
    NR=20;
    NP=20;
    NT=20;
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

    R_mu=(Rmax-Rmin)/2.;
    R_sd=R_mu/R_C;
    P_mu=(Pmax-Pmin)/2.;
    P_sd=Pmax/P_C;
    T_mu=(Tmax-Tmin)/2.;
    T_sd=Tmax/T_C;

    //for density
    T_mean = (Tmax-Tmin)/2.0; //in dimensionless units - this should be at a.
    T_stde = T_mean*0.35; //arbitrary
    P_mean = (Pmax-Pmin)/2.0; //in dimensionless units - this should be at a.
    P_stde = P_mean*0.35; //arbitrary
    R_mean = (Rmax-Rmin)/2.0; //in dimensionless units - this should be at a.
    R_stde = R_mean*0.35; //arbitrary

    //for counter thing
    mass = 1.;

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

    delete[] t_ana;
    t_ana = new double**[NT+1];
    for(int i =0; i<(NT+1); i++){
      t_ana[i] = new double*[NP+1];
      for(int j =0; j<(NP+1); j++){
        t_ana[i][j] = new double[NR+1];
        for(int k = 0; k<(NR+1); k++){
          t_ana[i][j][k] = 0.;
        }
      }
    }

    delete[] t_num;
    t_num = new double**[NT+1];
    for(int i =0; i<(NT+1); i++){
      t_num[i] = new double*[NP+1];
      for(int j =0; j<(NP+1); j++){
        t_num[i][j] = new double[NR+1];
        for(int k = 0; k<(NR+1); k++){
          t_num[i][j][k] = 0.;
        }
      }
    }

    delete[] total_mass;
    total_mass = new double**[NT+1];
    for(int i = 0; i<(NT+1); i++){
      total_mass[i] = new double*[NP+1];
      for(int j = 0; j<(NP+1); j++){
        total_mass[i][j] = new double[NR+1];
        for(int k = 0; k<(NR+1); k++){
          total_mass[i][j][k] = 0.;
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

    delete[] T_vec;
    T_vec = new double[noparticles];
    delete[] P_vec;
    P_vec = new double[noparticles];
    delete[] R_vec;
    R_vec = new double[noparticles];

    delete[] den;
    den = new double**[NT+1];
    for(int i=0; i<(NT+1); i++){
      den[i] = new double*[NP+1];
      for(int j=0; j<(NP+1); j++){
        den[i][j] = new double[NR+1];
        for(int k = 0; k<(NR+1); k++){
          den[i][j][k] = 0.;
        }
      }
    }

    //cout << "delete end" << endl;

    //adding initial conditions
    Ta[0] = 0.;
    dTb[0] = 1.;
    Pa[0] = 0.; //segmentation fault here ??
    dPb[0] = 1.;
    Ra[0] = 0.;
    //dRb[0] = 1.;
    //t[0] = 0.;

    //running code
    //cout << "starting" << endl;
    make_gauss(NT, NP, NR);
    //cout << "im working" << endl;
    find_inv(NT, NP, NR,T_B, P_B, R_B, T_g, P_g, R_g);
  //  cout << "still working" << endl;
    find_DR(NT, NP, NR, T_A=((Tmax-Tmin)/T_suminv), P_A=((Pmax-Pmin)/P_suminv), R_A=((Rmax-Rmin)/R_suminv), T_inv, P_inv, R_inv);
    //cout << "found that DR" << endl;
    // file_creator_gaussian(g,x);
    // file_creator_DR(DR,x);
    build_grid(NT, NP, NR, DT, DP, DR);
    //cout << "built ya grid" << endl;
    d_fill(NT, NP, NR,T_mean, P_mean, R_mean, T_stde, P_stde, R_stde);
    calculate_mass(NT,NP,NR,d);
    density_fill(NT, NP, NR,T_mean, P_mean, R_mean, T_stde, P_stde, R_stde);
    //cout << "filllliiiinggggggggg density" << endl;
    opacity_fill(NT, NP, NR);
    //cout << "okay filled opacity" << endl;
    calculate_optical_depth_ana(NT,NP,NR,kappa,d);
    //cout << "CALCULATED OPTICAL DEPTH YAS" << endl;
    calculate_optical_depth_num(NT,NP,NR,kappa,den);
    //cout << "CALCULATED OPTICAL DEPTH AGAIN YAS" << endl;

    //file creators
    file_creator_t_ana(NT,NP,NR,t_ana);
    file_creator_t_num(NT,NP,NR,t_num);
    //cout << "created file t" << endl;
    file_creator_phi(Pa);
    //cout << "created file phi" << endl;
    file_creator_theta(Ta);
    //cout << "created file theta" << endl;
    file_creator_d(d,NT,NP,NR);
    file_creator_den(den,NT,NP,NR);
    file_creator_T_vec(T_vec,NT,NP,NR);
    file_creator_P_vec(P_vec,NT,NP,NR);
    file_creator_R_vec(R_vec,NT,NP,NR);
    // //file_creator_gauss(gauss,Rb,NR);
    // file_creator_DR(DR,Ra);
    //file_creator_xRa(x,Ra);
    //file_creator_d(d,Ta,Pa,Ra,NT,NP,NR);
    //file_creator_den(den,Ta,Pa,Ra,NT,NP,NR);
    // file_creator_total_mass(NT,NP,NR,total_mass);
    //cout << "done" << endl;
    //cout << "BAA5" << endl;
  }

  return 0;
}
