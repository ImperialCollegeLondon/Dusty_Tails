#include "CAV.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include<numeric>

using namespace std;

//R GRID
//cell width
vector <double> make_gauss_R(double nr){
  for(double i=0.; i<=(C_R*sd_R); i+=(1./nr)){
    x_R.push_back(i);
    g_new_R = exp(-pow(i-mu_R,2.)/(2.*pow(sd_R,2.)));
    g_R.push_back(g_new_R);
  }
  return g_R;
  return x_R;
}

vector <double> find_inv_R(double B_R,vector<double> g_R){
  for(double i=0.; i<g_R.size(); i++){
    inv_new_R = 1.-(B_R*g_R[i]);
    if (inv_new_R == 0.){
      cout << "0 in DR vector" << endl; //so that DR isn't 0
    }
    inv_R.push_back(inv_new_R);
  }
  suminv_R = accumulate(inv_R.begin(),inv_R.end(),0.); //calculates sum of inverse gaussian values
  return inv_R;
}

vector <double> find_DR(double A_R, vector<double> inv_R){
  for(double i=0.; i<inv_R.size(); i++){
    g_new_R = A_R*inv_R[i];
    DR.push_back(g_new_R);
  }
  sumDR = accumulate(DR.begin(),DR.end(),0.);
  return DR;
}

int is(double NR){
  return 0;
}

int ie(double NR){
  return int(NR);
}

void build_grid_R(double NR, vector<double> DR){
  for (i=is(NR)+1; i<=ie(NR); i++){ // defines the cell edges
    Ra_new = Ra[i-1]+DR[i-1];
    Ra.push_back(Ra_new);
  }

  for (i=is(NR); i<=ie(NR)-1; i++){ // defines the cell centers
    Rb_new = Ra[i]+DR[i]/2.;
    Rb.push_back(Rb_new);
  }

  for (i=is(NR);i<=ie(NR)-1;i++){ // defines the width of a cell
    dRa_new = Ra[i+1]-Ra[i];
    dRa.push_back(dRa_new);
  }

  for (i=is(NR)+1;i<=ie(NR)-1;i++){ // defines the width between two cell centers
    dRb_new = Rb[i]-Rb[i-1];
    dRb.push_back(dRb_new);
  }
}



//THETA GRID
vector <double> make_gauss_T(double nt){
  for(double i=0.; i<=(C_T*sd_T); i+=(1./nt)){
    x_T.push_back(i);
    g_new_T = exp(-pow(i-mu_T,2.)/(2.*pow(sd_T,2.)));
    g_T.push_back(g_new_T);
  }
  return g_T;
  return x_T;
}

vector <double> find_inv_T(double B_T,vector<double> g_T){
  for(double i=0.; i<g_T.size(); i++){
    inv_new_T = 1.-(B_T*g_T[i]);
    if (inv_new_T == 0.){
      cout << "0 in DT vector" << endl; //so that DT isn't 0
    }
    inv_T.push_back(inv_new_T);
  }
  suminv_T = accumulate(inv_T.begin(),inv_T.end(),0.); //calculates sum of inverse gaussian values
  return inv_T;
}

vector <double> find_DT(double A_T, vector<double> inv_T){
  for(double i=0.; i<inv_T.size(); i++){
    g_new_T = A_T*inv_T[i];
    DT.push_back(g_new_T);
  }
  sumDT = accumulate(DT.begin(),DT.end(),0.);
  return DT;
}

int js(double NT){
  return 0;
}

int je(double NT){
  return int(NT);
}

void build_grid_T(double NT, vector<double> DT){
  for (i=js(NT)+1; i<=je(NT); i++){ // defines the cell edges
    Ta_new = Ta[i-1]+DT[i-1];
    Ta.push_back(Ta_new);
  }

  for (i=js(NT); i<=je(NT)-1; i++){ // defines the cell centers
    Tb_new = Ta[i]+DT[i]/2.;
    Tb.push_back(Tb_new);
  }

  for (i=js(NT);i<=je(NT)-1;i++){ // defines the width of a cell
    dTa_new = Ta[i+1]-Ta[i];
    dTa.push_back(dTa_new);
  }

  for (i=js(NT)+1;i<=je(NT)-1;i++){ // defines the width between two cell centers
    dTb_new = Tb[i]-Tb[i-1];
    dTb.push_back(dTb_new);
  }
}



//PHI GRID
vector <double> make_gauss_P(double np){
  for(double i=0.; i<=(C_P*sd_T); i+=(1./np)){
    x_P.push_back(i);
    g_new_P = exp(-pow(i-mu_P,2.)/(2.*pow(sd_P,2.)));
    g_P.push_back(g_new_P);
  }
  return g_P;
  return x_P;
}

vector <double> find_inv_P(double B_P,vector<double> g_P){
  for(double i=0.; i<g_P.size(); i++){
    inv_new_P = 1.-(B_P*g_P[i]);
    if (inv_new_P == 0.){
      cout << "0 in DP vector" << endl; //so that DT isn't 0
    }
    inv_P.push_back(inv_new_P);
  }
  suminv_P = accumulate(inv_P.begin(),inv_P.end(),0.); //calculates sum of inverse gaussian values
  return inv_P;
}

vector <double> find_DP(double A_P, vector<double> inv_P){
  for(double i=0.; i<inv_P.size(); i++){
    g_new_P = A_P*inv_P[i];
    DP.push_back(g_new_P);
  }
  sumDT = accumulate(DP.begin(),DP.end(),0.);
  return DP;
}

int ks(double NP){
  return 0;
}

int ke(double NP){
  return int(NP);
}

void build_grid_P(double NP, vector<double> DP){
  for (i=ks(NP)+1; i<=ke(NP); i++){ // defines the cell edges
    Pa_new = Pa[i-1]+DP[i-1];
    Pa.push_back(Pa_new);
  }

  for (i=ks(NP); i<=ke(NP)-1; i++){ // defines the cell centers
    Pb_new = Pa[i]+DP[i]/2.;
    Pb.push_back(Pb_new);
  }

  for (i=ks(NP);i<=ke(NP)-1;i++){ // defines the width of a cell
    dPa_new = Pa[i+1]-Pa[i];
    dPa.push_back(dPa_new);
  }

  for (i=ks(NP)+1;i<=ke(NP)-1;i++){ // defines the width between two cell centers
    dPb_new = Pb[i]-Pb[i-1];
    dPb.push_back(dPb_new);
  }
}



vector < vector < vector <double> > > d(NR, vector < vector <double> > (NT, vector<double>(NP)));

//fill density array
/*void density_fill(double NR, double NT, double NP)
{
  mean = 1.0; //mean is(NR) at a, which is(NR) 1 in dimensionless units
  stde = 0.2*mean;*/ //arbitrary

// Constructing density array
vector < vector < vector < double > > > density_fill_new(double NR, double NT, double NP){

  mean = 1.0; //mean is(NR) at a, which is(NR) 1 in dimensionless units
  stde = 0.2*mean;
  for(i=is(NR);i<=ie(NR)-1;i++){
    gauss_new_R = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Rb[i])-mean,2.)/(2.*pow(stde,2.)));
    gauss_R.push_back(gauss_new_R);
    cout << "BEEE" << endl;
    for(j=js(NT);j<=je(NT)-1;j++){
      gauss_new_T = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Tb[j])-mean,2.)/(2.*pow(stde,2.)));
      gauss_T.push_back(gauss_new_T);

      for(k=ks(NP);k<=ke(NP)-1;k++){
        gauss_new_P = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Pb[k])-mean,2.)/(2.*pow(stde,2.)));
        gauss_P.push_back(gauss_new_P);
        // cout << "BAAA2" << endl;

        vector < vector <double> > i_vector = d.at(i);
        vector<double> j_vector = i_vector.at(j);
        j_vector[k] = gauss_R[i]*gauss_T[j]*gauss_P[k];
        // d[i][j][k] = gauss_R[i]*gauss_T[j]*gauss_P[k];
      }
    }
  }
  return d;
}

//fill opacity array
/*void opacity_fill(double NR)
{
  //constructing opacity array
  for(i=is(NR); i<=ie(NR)-1; i++){
    k_new_R = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k_R.push_back(k_new_R);
  }

  for(i=js(NT); i<=je(NT)-1; i++){
    k_new_T = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k_T.push_back(k_new_T);
  }

  for(i=ks(NP); i<=ke(NP)-1; i++){
    k_new_P = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k_P.push_back(k_new_P);
  }

}

//calculate optical depth and fill array
void calculate_optical_depth(double NR)
{
  //constructing optical depth array
  for(i=is(NR)+1; i<=ie(NR)-1; i++){
    t_new = t[i-1]+(k[i]*d[i]*dRa[i]);
    t.push_back(t_new);
  }
  //cout << "done calculate optical depth" << endl;
}

//file creators
void file_creator_gaussian(vector<double> g, vector <double> x) {//textfile of gaussian vs Rb
  ofstream myfile ("gaussian_grid.txt");
  if (myfile.is_open()) {
    for (i=0.; i<g.size(); i++) {
      char string[15];

      myfile << x[i] << ",";
      myfile << g[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_inv(vector<double> DR, vector <double> x) {//textfile of non uniform grid DR vs Rb
  ofstream myfile ("inverse_gaussian_grid.txt");
  if (myfile.is_open()) {
    for (i=0.; i<DR.size(); i++) {
      char string[15];

      myfile << x[i] << ",";
      myfile << DR[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_DR(vector<double> DR, vector <double> Ra) {//textfile of density vs Rb
  ofstream myfile ("DR.txt");
  if (myfile.is_open()) {
    for (i=0.; i<DR.size(); i++) {
      char string[15];

      myfile << Ra[i] << ",";
      myfile << DR[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_xRa(vector<double> x, vector <double> Ra) {//textfile of density vs Rb
  ofstream myfile ("xRa.txt");
  if (myfile.is_open()) {
    for (i=0.; i<x.size(); i++) {
      char string[15];

      myfile << x[i] << ",";
      myfile << Ra[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_t(vector<double> t, vector <double> Ra, double NR) {//textfile of optical depth vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<t.size(); i++) {
      char string[15];

      myfile << Ra[i] << ",";
      myfile << t[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_gauss(vector<double> gauss, vector <double> Ra, double NR) {//textfile of density vs Rb
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

}
*/
