#include "file_creators_two.cpp"


int main() {
  semimajor(Period_days);
  //cout << a << endl;

  //parameters controlling the grid
  NR=200;
  NP=200;
  NT=200;

  rmin = 0.98*a;
  rmax = 1.13*a;
  Rmin = rmin/a;
  Rmax = rmax/a;

  Pmin = -0.7;
  Pmax = 0.01;

  Tmin = 1.55;
  Tmax = 1.65;

  R_B=0.2;
  P_B=0.2;
  T_B=0.2;

  // R_C=3.;
  // P_C=3.;
  // T_C=3.;

  R_mu=(Rmax-Rmin)/2.;
  R_sd=R_mu*0.1;
  P_mu=(Pmax-Pmin)/2.;
  P_sd=Pmax*0.1;
  T_mu=(Tmax-Tmin)/2.;
  T_sd=Tmax*0.1;

  //clearing pointers/vectors
  delete[] T_x;
  T_x = new double[NT+1];
  delete[] P_x;
  P_x = new double[NP+1];
  delete[] R_x;
  R_x = new double[NR+1];
  delete[] T_g;
  T_g = new double[NT+1];
  delete[] P_g;
  P_g = new double[NP+1];
  delete[] R_g;
  R_g = new double[NR+1];
  delete[] T_inv;
  T_inv = new double[NT+1];
  delete[] P_inv;
  P_inv = new double[NP+1];
  delete[] R_inv;
  R_inv = new double[NR+1];
  delete[] DT;
  DT = new double[NT+1];
  delete[] DP;
  DP = new double[NP+1];
  delete[] DR;
  DR = new double[NR+1];
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

  //adding initial conditions
  Ta[0] = Tmin;
  Pa[0] = Pmin;
  Ra[0] = Rmin;

  //building grid
  make_gauss(NT, NP, NR);
  find_inv(NT, NP, NR,T_B, P_B, R_B, T_g, P_g, R_g);
  find_DR(NT, NP, NR, T_A=((Tmax-Tmin)/T_suminv), P_A=((Pmax-Pmin)/P_suminv), R_A=((Rmax-Rmin)/R_suminv), T_inv, P_inv, R_inv);
  build_grid(NT, NP, NR, DT, DP, DR);

  //file creators
  // file_creator_phi(Pa);
  // file_creator_theta(Ta);
  // file_creator_R(Ra);

  //parameters controlling TIME
  TOLERANCE = 0.005;
  //looping over different orbital times to find optical depth
  for(TIME = 2.; TIME <= 2.; TIME=TIME+0.1){

    //clearing vectors and pointers
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

    data.reset();
    dataPointsByTime.clear();
    Particle.clear();
    points.clear();
    xPositions.clear();
    yPositions.clear();
    zPositions.clear();
    rPositions.clear();
    pPositions.clear();
    tPositions.clear();
    idPoints.clear();
    timePoints.clear();
    Mass.clear();
    xPositionsAtTime.clear();
    yPositionsAtTime.clear();
    zPositionsAtTime.clear();
    rPositionsAtTime.clear();
    pPositionsAtTime.clear();
    tPositionsAtTime.clear();
    idPointsAtTime.clear();
    MassAtTime.clear();
    SizeAtTime.clear();

    //filling grid and calculating optical depth
    get_positions("output_035.bin");
    density_fill(NT, NP, NR, tPositionsAtTime, pPositionsAtTime, rPositionsAtTime, MassAtTime, SizeAtTime, Ta, Pa, Ra);
    calculate_optical_depth(NT,NP,NR,kappa,den);

    //file creators
    // file_creator_rPositions(rPositionsAtTime);
    // file_creator_pPositions(pPositionsAtTime);
    // file_creator_tPositions(tPositionsAtTime);
    // file_creator_t(NT,NP,NR,t,TIME);
    // file_creator_den(den,NT,NP,NR,TIME);

    // cout << TIME << endl;
  }

  return 0;
}
