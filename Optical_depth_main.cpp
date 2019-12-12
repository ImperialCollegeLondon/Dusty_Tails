#include "Uniform_grid.cpp"


int main() {
  semimajor(Period_days);
  cout << "Semi-major axis = " << a << endl;

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

  for(NR=10.;NR<=1000.;(NR=NR+1.)){ //for uniform grid, calculating error
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
  }

  return 0;
}
