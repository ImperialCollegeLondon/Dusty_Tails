#include "Non_uniform_grid_1D_pointers.cpp"


int main() {
  semimajor(Period_days);
  cout << "Semi-major axis = " << a << endl;

  //setting constants
  NR=10;
  nr=(NR-1.)/2.;
  B=0.9; //A and B are in terms of a
  C=6.; //sd of gaussian grid = 2/C
  mu=1.;
  sd=Rmax/C;

  //clearing vectors
  delete[] Ra;
  Ra = new double[11];
  delete[] Rb;
  Rb = new double[11];
  delete[] dRa;
  dRa = new double[11];
  delete[] dRb;
  dRb = new double[11];
  delete[] d;
  d = new double[11];
  delete[] k;
  k = new double[11];
  delete[] t;
  t = new double[11];
  //gauss.clear();
  // g.clear();
  delete[] g;
  g = new double[11];
  delete[] DR;
  DR = new double[11];
  delete[] x;
  x = new double[11];
  delete[] inv;
  inv = new double[11];

  //adding initial conditions
  Ra[0] = 0.;
  dRb[0] = 1.;
  t[0] = 0.;

  //running code
  //cout << "starting" << endl;
  //cout << sd << endl;
  make_gauss();
  //cout << "done one" << endl;
  find_inv(B=0.9,g);
  find_DR(A=(2./suminv),inv);
  // file_creator_gaussian(g,x);
  // file_creator_DR(DR,x);
  build_grid(NR,DR);
  density_fill(NR);
  opacity_fill(NR);
  calculate_optical_depth(NR);

  //file creators
  // file_creator_t(t,Rb,NR);
  // //file_creator_gauss(gauss,Rb,NR);
  // file_creator_DR(DR,Ra);
  //file_creator_xRa(x,Ra);

  return 0;
}
