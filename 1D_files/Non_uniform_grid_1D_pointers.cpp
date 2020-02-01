#include "Constants.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include<numeric>

using namespace std;

//cell width
// vector <double> make_gauss(double nr){
//   for(double i=0.; i<=(C*sd); i+=(1./nr)){
//     x.push_back(i);
//     gauss_new = exp(-pow(i-mu,2.)/(2.*pow(sd,2.)));
//     g.push_back(gauss_new);
//     //cout << i << " " << gauss_new << endl; / prints the gaussian vector elements
//   }
//   //cout << "done make gauss" << endl;
//   return g;
//   // return x;
// }

double* make_gauss() {
  for(int i=0.; i<=10; i++){
    double interval = 1/5 * i;
    x[i] = interval;
    gauss_new = exp(-pow(interval-mu,2.)/(2.*pow(sd,2.)));
    g[i] = gauss_new;
  }
  return g;
}

double* find_inv(double B, double* g){
  double* inv = new double[11];
  double inv_sum = 0.;
  for(int i=0.; i<11; i++){
    inv_new = 1.-(B * g[i]);
    if (inv_new == 0.){
      cout << "0 in DR vector" << endl; //so that DR isn't 0
    }
    inv[i] = inv_new;
    inv_sum += inv_new;
  }
  suminv = inv_sum;
  //cout << suminv << endl; //prints the sum of the inverse gaussian values
  return inv;
}

double* find_DR(double A, double* inv) {
  double DR_sum = 0.;
  for(int i=0; i<11; i++){
    gauss_new = A*inv[i];
    DR[i] = gauss_new;
    DR_sum += gauss_new;
  }
  sumDR = DR_sum;
  //cout << sumDR << endl; //prints the sum of the DR values
  return DR;
}

//starting indices
int is(int NR){
  return 0;
}

//ending indices
int ie(int NR){
  return int(NR);
}

//building the grid
void build_grid(int NR, double* DR)
{
  for (i=is(NR)+1; i<=ie(NR); i++){ // defines the cell edges
    //cout << DR[i] << endl;
    Ra_new = Ra[i-1]+DR[i-1];
    Ra[i] = Ra_new;
    //cout << Ra[i] << endl;
  }

  for (i=is(NR); i<=ie(NR)-1; i++){ // defines the cell centers
    Rb_new = Ra[i]+DR[i]/2.;
    Rb[i] = Rb_new;
    //cout << Rb[i] << endl;
  }

  for (i=is(NR);i<=ie(NR)-1;i++){ // defines the width of a cell
    dRa_new = Ra[i+1]-Ra[i];
    dRa[i] = dRa_new;
    //cout << i << " " << dRa[i] << endl;
  }

  for (i=is(NR)+1;i<=ie(NR)-1;i++){ // defines the width between two cell centers
    dRb_new = Rb[i]-Rb[i-1];
    dRb[i] = dRb_new;
    //cout << dRb[i] << endl;
  }
}

//fill density array
void density_fill(int NR)
{
  mean = 1.0; //mean is(NR) at a, which is(NR) 1 in dimensionless units
  stde = 0.2*mean; //arbitrary

  // Constructing density array
  for(i=is(NR); i<=ie(NR)-1; i++){
    d_new = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Rb[i])-mean,2.)/(2.*pow(stde,2.)));
    d[i] = d_new;
  }
}

//fill opacity array
void opacity_fill(int NR)
{
  //constructing opacity array
  for(i=is(NR); i<=ie(NR)-1; i++){
    k_new = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k[i] = k_new;
  }
  //cout << "done opacity fill" << endl;
}

//calculate optical depth and fill array
void calculate_optical_depth(int NR)
{
  //constructing optical depth array
  for(i=is(NR)+1; i<=ie(NR)-1; i++){
    t_new = t[i-1]+(k[i]*d[i]*dRa[i]);
    t[i] = t_new;
  }
  //cout << "done calculate optical depth" << endl;
}

//file creators
// void file_creator_gaussian(vector<double> g, vector <double> x) {//textfile of gaussian vs Rb
//   ofstream myfile ("gaussian_grid.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<g.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << g[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }
//
// void file_creator_inv(vector<double> DR, vector <double> x) {//textfile of non uniform grid DR vs Rb
//   ofstream myfile ("inverse_gaussian_grid.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<DR.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << DR[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }
//
// void file_creator_DR(vector<double> DR, vector <double> Ra) {//textfile of density vs Rb
//   ofstream myfile ("DR.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<DR.size(); i++) {
//       char string[15];
//
//       myfile << Ra[i] << ",";
//       myfile << DR[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

// void file_creator_xRa(vector<double> x, vector <double> Ra) {//textfile of density vs Rb
//   ofstream myfile ("xRa.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<x.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << Ra[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

void file_creator_t(double* t, double* Ra, int NR) {//textfile of optical depth vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<11; i++) {
      char string[15];

      myfile << Ra[i] << ",";
      myfile << t[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

/*void file_creator_gauss(vector<double> gauss, vector <double> Ra, int NR) {//textfile of density vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<gauss.size(); i++) {
      char string[15];

      myfile << Ra[i] << ",";
      myfile << gauss[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}*/
