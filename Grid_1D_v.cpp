#include "Grid_1D_v.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include "constants.h"

using namespace std;

void build_grid(void)
{
  int       i;
  double    DR;

  // Constructing 1D row of cells with constant width, DR (for uniform grid only)
  DR = (Rmax-Rmin) / NR;

  // Asign starting indicies
  is = 0;
  ie = NR;
  itot = NR + 1; //number of a-mesh indices

  // Buldng the grid
  for (i=is+1; i<=ie; i++){ // defines the cell edges
    Ra_new = Ra[i-1]+DR;
    Ra.push_back(Ra_new);
  }

  for (i=is; i<=ie-1; i++){ // defines the cell centers
    Rb_new = Ra[i]+DR/2.;
    Rb.push_back(Rb_new);
  }

  for (i=is;i<=ie-1;i++){ // defines the width of a cell
    dRa_new = Ra[i+1]-Ra[i];
    dRa.push_back(dRa_new);
  }

  for (i=is+1;i<=ie-1;i++){ // defines the width between two cell centers
    dRb_new = Rb[i]-Rb[i-1];
    dRb.push_back(dRb_new);
  }


cout << DR << endl;

}

void density_fill(void) //fill density vector
{ int i;

  //boundaries
  is=0;
  ie=NR;
  itot=NR+1;

  // Constructing density array
  for(i=is; i<=ie-1; i++){
    d.push_back(0.85);
  }

}

void opacity_fill(void) //fill opacity vector
{ int i;

  //boundaries
  is=0;
  ie=NR;
  itot=NR+1;

  //constructing opacity array
  for(i=is; i<=ie-1; i++){
    k.push_back(6.1e14);
  }

}

void calculate_optical_depth()
{ int i;

  //boundaries
  is=0;
  ie=NR;
  itot=NR+1;

  //constructing optical depth array
  for(i=is+1; i<=ie-1; i++){
    t_new = t[i-1]+(k[i]*d[i]*dRa[i]);
    t.push_back(t_new);
    cout << t[i] << endl;
  }

}

void file_creator(vector<double> t, vector<double> Rb) {
  ofstream myfile ("optical_depth.txt");
  if (myfile.is_open()) {
    for (int i=0.; i<t.size(); i++) {
      char string[15];

      myfile << Rb[i] << ",";
      myfile << t[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

int main() {

  Ra.push_back(0.);
  dRb.push_back(1.);
  build_grid();

  density_fill();

  opacity_fill();

  t.push_back(0.);
  calculate_optical_depth();

  file_creator(t,Rb);

  return 0;
}
