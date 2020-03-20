#include "Non_uniform_grid_3D.cpp"


int main() {
  semimajor(Period_days);
  cout << "Semi-major axis = " << a << endl;

  //for non uniform 3D grid
  //constants
  NR = 10.;
  NT = 10.;
  NP = 10.;
  nr=(NR-1.)/2.;
  nt=(NT-1.)/2.;
  np=(NP-1.)/2.;
  B_R=0.9; //A and B are in terms of a
  B_T=0.9;
  B_P=0.9;

  C_R=6.; //sd of gaussian grid = 2/C
  C_T=6.;
  C_P=6.;

  mu_R=1.;
  sd_R=Rmax/C_R;
  mu_T=1.;
  sd_R=Tmax/C_T;
  mu_P=1.;
  sd_T=Pmax/C_P;

  //clearing vectors
  Ra.clear();
  Rb.clear();
  dRb.clear();
  dRa.clear();
  Ta.clear();
  Tb.clear();
  dTa.clear();
  dTb.clear();
  Pa.clear();
  Pb.clear();
  dPa.clear();
  dPb.clear();
  gauss_R.clear();
  gauss_T.clear();
  gauss_P.clear();
  g_R.clear();
  g_T.clear();
  g_P.clear();
  x_R.clear();
  x_T.clear();
  x_P.clear();
  d.clear();
  //k.clear();
  t.clear();
  DR.clear();
  DP.clear();
  DT.clear();
  inv_R.clear();
  inv_T.clear();
  inv_P.clear();

  //adding initial conditions
  Ra.push_back(0.);
  dRb.push_back(1.);
  t.push_back(0.);

  //running code
  make_gauss_R(nr);
  find_inv_R(B_R,g_R);
  find_DR(A_R=(2./suminv_R),inv_R);
  build_grid_R(NR,DR);
  cout << "BAA1" << endl;

  make_gauss_T(nt);
  find_inv_T(B_T,g_T);
  find_DT(A_T=(2./suminv_T),inv_T);
  build_grid_T(NT,DT);
  cout << "BAA2" << endl;

  make_gauss_P(np);
  find_inv_P(B_P,g_P);
  find_DP(A_P=(2./suminv_P),inv_P);
  build_grid_P(NP,DP);
  cout << "BAA3" << endl;

  density_fill_new(NR,NT,NP);
  cout << d[i][0][0] << "," << d[0][j][0] << "," << d[0][0][k] << endl;
  //opacity_fill(NR);
  //calculate_optical_depth(NR);

  //file creators
  //file_creator_t(t,Rb,NR);
  //file_creator_gauss(gauss,Rb,NR);
  //file_creator_DR(DR,Ra);
  //file_creator_xRa(x,Ra);

  return 0;
}
