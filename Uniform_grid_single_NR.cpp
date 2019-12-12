#include "CAV.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>

using namespace std;

double DR;
int is=0.;
int ie=NR;


//building the grid
void build_grid(void)
{
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
}

//fill density array
void density_fill(void)
{
  mean = 1.0; //mean is at a, which is 1 in dimensionless units
  sd = 0.2*mean; //arbitrary

  // Constructing density array
  for(i=is; i<=ie-1; i++){
    gauss_new = (density*pow(a,3.)/Mstar_kg)*(1./(sd*sqrt(2*M_PI))*exp(-pow((Rb[i])-mean,2.)/(2.*pow(sd,2.))));
    gauss.push_back(gauss_new);
    d.push_back(gauss_new);
  }
}

//fill opacity array
void opacity_fill(void)
{
  //constructing opacity array
  for(i=is; i<=ie-1; i++){
    k_new = (3./4.)*(1./density)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k.push_back(k_new);
  }
}

//calculate optical depth and fill array
void calculate_optical_depth()
{
  //constructing optical depth array
  for(i=is+1; i<=ie-1; i++){
    t_new = t[i-1]+(k[i]*d[i]*dRa[i]);
    t.push_back(t_new);
  }
}

//file creators
void file_creator_t(vector<double> t, vector <double> Rb) {//textfile of optical depth vs Rb
  ofstream myfile ("optical_depth.txt");
  if (myfile.is_open()) {
    for (i=0.; i<t.size(); i++) {
      char string[15];

      myfile << Rb[i] << ",";
      myfile << t[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_gauss(vector<double> gauss, vector <double> Rb) {//textfile of density vs Rb
  ofstream myfile ("gaussian.txt");
  if (myfile.is_open()) {
    for (i=0.; i<gauss.size(); i++) {
      char string[15];

      myfile << Rb[i] << ",";
      myfile << gauss[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}



/*void looping(void){
  for(NR=10.;NR<=1000.;NR=NR*10.){
    build_grid();
    density_fill();
    opacity_fill();
    calculate_optical_depth();
    file_creator_gauss(gauss,Rb);
    file_creator_t(t,Rb);
  }
}*/
