#include "3D_files/Const.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <array>
#include <chrono>

using namespace std;

double* vec;
double* den;

int NR = 50;
double Rmax = 10.0;
double Rmin = 0.0;

double DR(double NR){
  return (Rmax-Rmin)/NR;
}

//starting indices
int is(double NR){
  return 0;
}

//ending indices
int ie(double NR){
  return int(NR);
}

void build_grid(double NR)
{
  for (i=is(NR)+1; i<=ie(NR); i++){ // defines the cell edges
    Ra_new = Ra[i-1]+DR(NR);
    Ra[i] = Ra_new;
  }

  for (i=is(NR); i<=ie(NR)-1; i++){ // defines the cell centers
    Rb_new = Ra[i]+DR(NR)/2.;
    Rb[i] = Rb_new;
  }

  for (i=is(NR);i<=ie(NR)-1;i++){ // defines the width of a cell
    dRa_new = Ra[i+1]-Ra[i];
    dRa[i] = dRa_new;
  }

  for (i=is(NR)+1;i<=ie(NR)-1;i++){ // defines the width between two cell centers
    dRb_new = Rb[i]-Rb[i-1];
    dRb[i] = dRb_new;
  }
}

int main(){
  cout << "BAA" << endl;
  delete[] Ra;
  Ra = new double[NR+1];
  delete[] Rb;
  Rb = new double[NR+1];
  delete[] dRa;
  dRa = new double[NR+1];
  delete[] dRb;
  dRb = new double[NR+1];


  Ra[0] = 0.;
  dRb[0] = 1.;


  DR(NR);
  build_grid(NR);

  const int mean = int((Rmax-Rmin)/2.0);
  const int stde = int(mean*0.35);

  double mass = 1.; //mass of particles
  int noparticles = 100000;

  // construct a trivial random generator engine from a time-based seed:
  //unsigned seedR = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seedR = 1;
  default_random_engine generatorR (seedR);
  normal_distribution<double> distributionR (mean,stde);


  //double arr[] = {5.5,4.5,5.5,5.5,6.5,5.5,8.5,2.5,5.5,4.5};
  //vector<double> vec(arr, arr + sizeof(arr) / sizeof(arr[0]) );

  delete[] vec;
  vec = new double[noparticles];
  delete[] den;
  den = new double[NR+1];


  for (int i=0; i<noparticles; i++){
    double vec_new = distributionR(generatorR);
    vec[i] = vec_new;
    //cout << vec[i] << endl;
  }

  for(int j=0;j<=NR;j++){
    double den_new = 0.;
    den[j] = den_new; //setting initial density to 0 in grid
  }

  for(int i=0; i<noparticles; i++){//iterating through particles
    for(int j=0; j<=NR; j++){//iterating through grid cells
      if((Ra[j]<=vec[i]) && (vec[i]<Ra[j+1])){
        den[j] = (den[j]+1);
      }
      else{}
    }
  }

  for(j=0;j<=NR;j++){
    cout << Ra[j] << endl;
  }

  for(j=0;j<=NR;j++){
    cout << den[j] << endl;
  }

  return 0;

}


// int j,k;
// vector < vector < vector < double > > > d;
//
// int is(double NR){
//   return 0;
// }
//
// int ie(double NR){
//   return int(NR);
// }
//
// int js(double NT){
//   return 0;
// }
//
// int je(double NT){
//   return int(NT);
// }
//
// int ks(double NP){
//   return 0;
// }
//
// int ke(double NP){
//   return int(NP);
// }
//
// vector<vector<vector<double>>> density_fill_new(double NR, double NT, double NP){
//   for(i=is(NR);i<=ie(NR)-1;i++){
//     gauss_new_R = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Rb[i])-mean,2.)/(2.*pow(stde,2.)));
//     gauss_R.push_back(gauss_new_R);
//     for(j=js(NT);j<=je(NT)-1;i++){
//       gauss_new_T = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Tb[i])-mean,2.)/(2.*pow(stde,2.)));
//       gauss_T.push_back(gauss_new_T);
//       for(k=ks(NP);k<=ke(NP)-1;i++){
//         gauss_new_P = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Pb[i])-mean,2.)/(2.*pow(stde,2.)));
//         gauss_P.push_back(gauss_new_P);
//         d[i][j][k] = gauss_R[i]*gauss_T[j]*gauss_P[k];
//         cout << d[i] << "," << d[j] << "," << d[k] << endl;
//       }
//     }
//   }
//   return d;
// }
//
// int main(){
//   density_fill_new(NR=10.,NT=10.,NP=10.);
//   return 0;
// }
