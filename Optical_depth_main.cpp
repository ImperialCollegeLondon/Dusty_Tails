#include "Non_uniform_grid.cpp"


int main() {
  semimajor(Period_days);
  cout << "Semi-major axis = " << a << endl;

  //for uniform grid, single NR, calculating error
  /*NR = 200.;
  Ra.push_back(0.);
  dRb.push_back(1.);
  t.push_back(0.);
  //find_DR(NR);
  build_grid(NR);
  density_fill(NR);
  opacity_fill(NR);
  calculate_optical_depth(NR);
  file_creator_t(t,Rb,NR);
  file_creator_gauss(gauss,Rb,NR);*/

  /*for(NR=10.;NR<=1000.;(NR=NR+1.)){ //for uniform grid, calculating error
    Ra.push_back(0.);
    dRb.push_back(1.);
    t.push_back(0.);
    //find_DR(NR);
    build_grid(NR);
    density_fill(NR);
    opacity_fill(NR);
    calculate_optical_depth(NR);
    file_creator_t(t,Rb,NR);
    file_creator_gauss(gauss,Rb,NR);
    Ra.clear();
    Rb.clear();
    dRb.clear();
    dRa.clear();
    d.clear();
    k.clear();
    t.clear();
    gauss.clear();
  }*/

  //for non-uniform grid
  NR=100.; //setting parameters
  A=0.03; //A and B are in terms of a
  B=0.001;
  Ra.clear();
  Rb.clear();
  dRb.clear();
  dRa.clear();
  d.clear();
  k.clear();
  t.clear();
  gauss.clear();
  g.clear();
  DR.clear();
  Ra.push_back(0.);
  dRb.push_back(1.);
  t.push_back(0.);
  make_gauss(NR);
  find_DR(A,B,g);
  file_creator_gaussian(g,x);
  file_creator_inv(DR,x);
  build_grid(NR,DR);
  density_fill(NR);
  opacity_fill(NR);
  calculate_optical_depth(NR);
  file_creator_t(t,Rb,NR);
  file_creator_gauss(gauss,Rb,NR);
  file_creator_DR(DR,Ra);
  file_creator_xRa(x,Ra);

  return 0;
}
