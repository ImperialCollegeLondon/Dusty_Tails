#include "Non_uniform_grid_3D.cpp"


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

  /*for(NR=11.;NR<=20.;(NR=NR+1.)){ //for uniform grid, multiple NR, calculating error
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

  //setting constants
  /*for(NR=11.;NR<=1000.;(NR=NR+1.)){
    nr=(NR-1.)/2.;
    B=0.9; //A and B are in terms of a
    C=6.; //sd of gaussian grid = 2/C
    mu=1.;
    sd=Rmax/C;

    //clearing vectors
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
    x.clear();
    inv.clear();

    //adding initial conditions
    Ra.push_back(0.);
    dRb.push_back(1.);
    t.push_back(0.);

    //running code
    //cout << "starting" << endl;
    //cout << sd << endl;
    make_gauss(nr);
    //cout << "done one" << endl;
    find_inv(B=0.9,g);
    find_DR(A=(2./suminv),inv);
    file_creator_gaussian(g,x);
    file_creator_DR(DR,x);
    build_grid(NR,DR);
    density_fill(NR);
    opacity_fill(NR);
    calculate_optical_depth(NR);

    //file creators
    file_creator_t(t,Rb,NR);
    file_creator_gauss(gauss,Rb,NR);
    file_creator_DR(DR,Ra);
    file_creator_xRa(x,Ra);
  }*/

  //for non uniform 3D grid
  //constants
  for(NR=11.;NR<=1000.;(NR=NR+1.)){
    nr=(NR-1.)/2.;
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
    d.clear();
    k.clear();
    t.clear();
    gauss.clear();
    g.clear();
    DR.clear();
    x.clear();
    inv.clear();

    //adding initial conditions
    Ra.push_back(0.);
    dRb.push_back(1.);
    t.push_back(0.);

    //running code
    //cout << "starting" << endl;
    //cout << sd << endl;
    make_gauss(nr);
    //cout << "done one" << endl;
    find_inv(B=0.9,g);
    find_DR(A=(2./suminv),inv);
    file_creator_gaussian(g,x);
    file_creator_DR(DR,x);
    build_grid(NR,DR);
    density_fill(NR);
    opacity_fill(NR);
    calculate_optical_depth(NR);

    //file creators
    file_creator_t(t,Rb,NR);
    file_creator_gauss(gauss,Rb,NR);
    file_creator_DR(DR,Ra);
    file_creator_xRa(x,Ra);
  }

  return 0;
}
