#include "Const.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include<numeric>

using namespace std;

//indices
int is(int NT){
  return 0;
}
int ie(int NT){
  return int(NT);
}

int js(int NP){
  return 0;
}
int je(int NP){
  return int(NP);
}

int ks(int NR){
  return 0;
}
int ke(int NR){
  return int(NR);
}

//cell width
void make_gauss(int NT, int NP, int NR) {
  for(int i=is(NT); i<=ie(NT); i++){
    double interval = 1/5 * i;
    T_x[i] = interval;
    double gauss_new = exp(-pow(interval-T_mu,2.)/(2.*pow(T_sd,2.)));
    T_g[i] = gauss_new;
  }
  for(j=js(NP); j<=je(NP); j++){
    double interval = 1/5 * j;
    P_x[j] = interval;
    double gauss_new = exp(-pow(interval-P_mu,2.)/(2.*pow(P_sd,2.)));
    P_g[j] = gauss_new;
  }
  for(k=ks(NR); k<=ke(NR); k++){
    double interval = 1/5 * k;
    R_x[k] = interval;
    double gauss_new = exp(-pow(interval-R_mu,2.)/(2.*pow(R_sd,2.)));
    R_g[k] = gauss_new;
  }
}


void find_inv(int NT, int NP, int NR, double T_B, double P_B, double R_B, double* T_g, double* P_g, double* R_g){
  double T_inv_sum = 0.;
  double P_inv_sum = 0.;
  double R_inv_sum = 0.;
  for(i=is(NT); i<=ie(NT); i++){
    double inv_new = 1.-(T_B * T_g[i]);
    if (inv_new == 0.){
      cout << "0 in DT vector" << endl; //so that DR isn't 0
    }
    T_inv[i] = inv_new;
    T_inv_sum += inv_new;
  }
  for(j=js(NP); j<=je(NP); j++){
    double inv_new = 1.-(P_B * P_g[j]);
    if (inv_new == 0.){
      cout << "0 in DP vector" << endl; //so that DR isn't 0
    }
    P_inv[j] = inv_new;
    P_inv_sum += inv_new;
  }
  for(k=ks(NR); k<=ke(NR); k++){
    double inv_new = 1.-(R_B * R_g[k]);
    if (inv_new == 0.){
      cout << "0 in DR vector" << endl; //so that DR isn't 0
    }
    R_inv[k] = inv_new;
    R_inv_sum += inv_new;
  }
  T_suminv = T_inv_sum;
  P_suminv = P_inv_sum;
  R_suminv = R_inv_sum;
}


void find_DR(int NT, int NP, int NR, double T_A, double P_A, double R_A, double* T_inv, double* P_inv, double* R_inv) {
  double DT_sum = 0.;
  double DP_sum = 0.;
  double DR_sum = 0.;
  for(i=is(NT); i<=ie(NT); i++){
    double gauss_new = T_A*T_inv[i];
    DT[i] = gauss_new;
    DT_sum += gauss_new;
  }
  for(j=js(NP); j<=je(NP); j++){
    double gauss_new = P_A*P_inv[j];
    DP[j] = gauss_new;
    DP_sum += gauss_new;
  }
  for(k=ks(NR); k<=ke(NR); k++){
    double gauss_new = R_A*R_inv[k];
    DR[k] = gauss_new;
    DR_sum += gauss_new;
  }
  sumDT = DT_sum;
  sumDP = DP_sum;
  sumDR = DR_sum;
}

//building the grid
void build_grid(int NT, int NP, int NR, double* DT, double* DP, double* DR){
//R grid
  for (i=is(NT)+1; i<=ie(NT); i++){ // defines the cell edges
    Ta_new = Ta[i-1]+DT[i-1];
    Ta[i] = Ta_new;
  }
  for (i=is(NT); i<=ie(NT)-1; i++){ // defines the cell centers
    Tb_new = Ta[i]+DT[i]/2.;
    Tb[i] = Tb_new;
  }
  for (i=is(NT);i<=ie(NT)-1;i++){ // defines the width of a cell
    dTa_new = Ta[i+1]-Ta[i];
    dTa[i] = dTa_new;
  }
  for (i=is(NT)+1;i<=ie(NT)-1;i++){ // defines the width between two cell centers
    dTb_new = Tb[i]-Tb[i-1];
    dTb[i] = dTb_new;
  }

//PHI grid
  for (j=js(NP)+1; j<=je(NP); j++){ // defines the cell edges
    Pa_new = Pa[j-1]+DP[j-1];
    Pa[j] = Pa_new;
  }
  for (j=js(NP); j<=je(NP)-1; j++){ // defines the cell centers
    Pb_new = Pa[j]+DP[j]/2.;
    Pb[j] = Pb_new;
  }
  for (j=js(NP);j<=je(NP)-1;j++){ // defines the width of a cell
    dPa_new = Pa[j+1]-Pa[j];
    dPa[j] = dPa_new;
  }
  for (j=js(NP)+1;j<=je(NP)-1;j++){ // defines the width between two cell centers
    dPb_new = Pb[j]-Pb[j-1];
    dPb[j] = dPb_new;
  }

//THETA grid
  for (k=ks(NR)+1; k<=ke(NR); k++){ // defines the cell edges
    Ra_new = Ra[k-1]+DR[k-1];
    Ra[k] = Ra_new;
  }
  for (k=ks(NR); k<=ke(NR)-1; k++){ // defines the cell centers
    Rb_new = Ra[k]+DR[k]/2.;
    Rb[k] = Rb_new;
  }
  for (k=ks(NR);k<=ke(NR)-1;k++){ // defines the width of a cell
    dRa_new = Ra[k+1]-Ra[k];
    dRa[k] = dRa_new;
  }
  for (k=ks(NR)+1;k<=ke(NR)-1;k++){ // defines the width between two cell centers
    dRb_new = Rb[k]-Rb[k-1];
    dRb[k] = dRb_new;
  }
}

//fill density array
void density_fill(int NT, int NP, int NR, float mean, float stde)
{
  //mean is(NR) at a, which is(NR) 1 in dimensionless units

  for(i=is(NT); i<=ie(NT); i++){
    for(j=js(NP); j<=je(NP); j++){
      for(k=ks(NR); k<=ke(NR); k++){
        d_new_T = exp(-pow((Tb[i])-mean,2.)/(2.*pow(stde,2.)));
        double* T_d = new double[51];
        T_d[i] = d_new_T;
        d_new_P = exp(-pow((Pb[j])-mean,2.)/(2.*pow(stde,2.)));
        double* P_d = new double[51];
        P_d[j] = d_new_P;
        d_new_R = exp(-pow((Rb[k])-mean,2.)/(2.*pow(stde,2.)));
        double* R_d = new double[51];
        R_d[k] = d_new_R;
        d[i][j][k] = (density*pow(a,3.)/Mstar_kg)*T_d[i]*P_d[j]*R_d[k];
      }
    }
  }
}

//fill opacity array
void opacity_fill(int NT, int NP, int NR)
{
  for(i=is(NT); i<=ie(NT); i++){
    for(j=js(NP); j<=je(NP); j++){
      for(k=ks(NR); k<=ke(NR); k++){
        kappa[i][j][k] = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
      }
    }
  }
  //cout << "done opacity fill" << endl;
}

//calculate optical depth and fill array
void calculate_optical_depth(int NT, int NP, int NR, double*** kappa, double*** d)
{
  //constructing optical depth array
  for(i=is(NT); i<=ie(NT); i++){
    for(j=js(NP); j<=je(NP); j++){
      for(k=ks(NR); k<=ke(NR); k++){
        t[i][j][k]=t[i][j][k-1]+(kappa[i][j][k]*d[i][j][k]*dRa[k]);
        cout << t[i][j][k] << endl;
      }
    }
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

void file_creator_t(double*** t) {//textfile of optical depth vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<51; i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for ( int j=0; j<51; j++){ //looping over PHI coordinates
        for (int k=0; k<51; k++){ //looping over R coordinates
          myfile << t[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_d(double*** d, double* Ta, double* Pa, double* Ra, int NT, int NP, int NR) {//textfile of density vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/density_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<51; i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for ( int j=0; j<51; j++){ //looping over PHI coordinates
        for (int k=0; k<51; k++){ //looping over R coordinates
          myfile << d[i][j][k] << endl; //d[theta][phi][R]
        }
      }
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
