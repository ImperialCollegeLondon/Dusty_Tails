#include "CAV.h"
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

using namespace std;

int main()
{
  const int nparticles = 100;

  double Rmax = 10.0;
  double Rmin = 0.0;

  int mean = 5.0;
  int stde = 2.0;

  default_random_engine generator;
  normal_distribution<double> gaussian(mean,stde);

  int p[10]={};

  for (int i=0; i<100; i++) { //100 experiments
    double number = gaussian(generator);
    if ((number>=Rmin)&&(number<Rmax)) ++p[int(number)]; //final bit increments in p and then adds the value
  }

  for(int i=0; i<10; i++){
    cout << p[i] << endl;
  }

  cout << "normal_distribution (5.0,2.0):" << endl;

  for (int i=0; i<10; i++) { //iterating through number of particles (10)
    //cout << p[i] << endl;
    //cout << i << "-" << (i+1) << ": ";
    //cout << string(p[i],'*') << endl;
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
