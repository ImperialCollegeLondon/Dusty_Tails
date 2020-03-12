#include "Const.h"
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
  for(int i=is(NT); i<=ie(NT-1); i++){
    double interval = (Tmax-Tmin)/(NT-1) * double(i);
    T_x[i] = interval;
    double gauss_new = exp(-pow(T_x[i]-T_mu,2.)/(2.*pow(T_sd,2.)));
    T_g[i] = gauss_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double interval = (Pmax-Pmin)/(NP-1) * double(j);
    P_x[j] = interval;
    double gauss_new = exp(-pow(P_x[j]-P_mu,2.)/(2.*pow(P_sd,2.)));
    P_g[j] = gauss_new;
  }
  for(int k=ks(NR); k<=ke(NR-1); k++){
    double interval = (Rmax-Rmin)/(NR-1) * double(k);
    R_x[k] = interval;
    double gauss_new = exp(-pow((R_x[k]-R_mu),2.)/(2.*pow(R_sd,2.)));
    R_g[k] = gauss_new;
  }
}


void find_inv(int NT, int NP, int NR, double T_B, double P_B, double R_B, double* T_g, double* P_g, double* R_g){
  double T_inv_sum = 0.;
  double P_inv_sum = 0.;
  double R_inv_sum = 0.;
  for(i=is(NT); i<=ie(NT-1); i++){
    double inv_new = 1.-(T_B * T_g[i]);
    if (inv_new == 0.){
      cout << "0 in DT vector" << endl; //so that DR isn't 0
    }
    T_inv[i] = inv_new;
    T_inv_sum += inv_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double inv_new = 1.-(P_B * P_g[j]);
    if (inv_new == 0.){
      cout << "0 in DP vector" << endl; //so that DR isn't 0
    }
    P_inv[j] = inv_new;
    P_inv_sum += inv_new;
  }
  for(k=ks(NR); k<=ke(NR-1); k++){
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
  for(i=is(NT); i<=ie(NT-1); i++){
    double gauss_new = T_A*T_inv[i];
    DT[i] = gauss_new;
    DT_sum += gauss_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double gauss_new = P_A*P_inv[j];
    DP[j] = gauss_new;
    DP_sum += gauss_new;
  }
  for(k=ks(NR); k<=ke(NR-1); k++){
    double gauss_new = R_A*R_inv[k];
    DR[k] = gauss_new;
    DR_sum += gauss_new;
    //cout << DR[k] << endl;
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

void d_fill(int NT, int NP, int NR, double T_mean, double P_mean, double R_mean, double T_stde, double P_stde, double R_stde)
{
  for(i=is(NT); i<=ie(NT); i++){
    for(j=js(NP); j<=je(NP); j++){
      for(k=ks(NR); k<=ke(NR); k++){
        d_new_T = exp(-pow((Tb[i])-T_mean,2.)/(2.*pow(T_stde,2.)));
        double* T_d = new double[NT+1];
        T_d[i] = d_new_T;
        d_new_P = exp(-pow((Pb[j])-P_mean,2.)/(2.*pow(P_stde,2.)));
        double* P_d = new double[NP+1];
        P_d[j] = d_new_P;
        d_new_R = exp(-pow((Rb[k])-R_mean,2.)/(2.*pow(R_stde,2.)));
        double* R_d = new double[NR+1];
        R_d[k] = d_new_R;
        d[i][j][k] = (density*pow(a,3.)/Mstar_kg)*T_d[i]*P_d[j]*R_d[k];
      }
    }
  }
}

double calculate_mass(int NT, int NP, int NR, double*** d){
  TMASS = 0.;
  for(i=is(NT); i<=ie(NT); i++){
    for(j=js(NP); j<=je(NP); j++){
      for(k=ks(NR); k<=ke(NR); k++){
        total_mass[i][j][k]=(d[i][j][k]*(((pow(Ra[k+1],3.)-pow(Ra[k],3.))/3.)*(-cos(Ta[i+1])+cos(Ta[i]))*(dPa[j])));
        double mass_new = total_mass[i][j][k]*(Mstar_kg);
        TMASS += mass_new;
      }
    }
  }
  cout << "TMASS = " << TMASS << endl;
  return TMASS;
}

void density_fill(int NT, int NP, int NR, double T_mean, double P_mean, double R_mean, double T_stde, double P_stde, double R_stde)
{
  unsigned seedT = 1;
  unsigned seedP = 2;
  unsigned seedR = 3;

  default_random_engine generatorT(seedT);
  normal_distribution<double> distributionT(T_mean,T_stde);

  default_random_engine generatorP(seedP);
  normal_distribution<double> distributionP(P_mean,P_stde);

  default_random_engine generatorR(seedR);
  normal_distribution<double> distributionR(R_mean, R_stde);

  for (int x=0; x<noparticles; x++){ //pointer which holds the particle positions in R
    double vec_new = distributionT(generatorT);
    T_vec[x] = vec_new;
    //cout << T_vec[x] << endl;
  }
  for (int x=0; x<noparticles; x++){ //pointer which holds the particle positions in R
    double vec_new = distributionP(generatorP);
    P_vec[x] = vec_new;
    //cout << P_vec[x] << endl;
  }
  for (int x=0; x<noparticles; x++){ //pointer which holds the particle positions in R
    double vec_new = distributionR(generatorR);
    R_vec[x] = vec_new;
    //cout << R_vec[x] << endl;
  }

  int p[20] = {0};
  int t[20] = {0};
  int r[20] = {0};
  tmass = 0.;

  for(int x=0; x<=noparticles; x++){//iterating through particles
    for(int i=is(NT)+1; i<ie(NT); i++){
      if((Ta[i]<=T_vec[x]) && (T_vec[x]<Ta[i+1])){
        t[i] += 1.;
        for(int j=js(NP)+1; j<je(NP); j++){
          if((Pa[j]<=P_vec[x]) && (P_vec[x]<Pa[j+1])){
            p[j] += 1.;
            for(int k=ks(NR)+1; k<ke(NR); k++){
              if((Ra[k]<=R_vec[x]) && (R_vec[x]<Ra[k+1])){
                r[k] += 1.;
                double den_new = den[i][j][k]+((TMASS/(Mstar_kg*noparticles))/(((pow(Ra[k+1],3.)-pow(Ra[k],3.))/3.)*(-cos(Ta[i+1])+cos(Ta[i]))*(dPa[j])));
                den[i][j][k] = den_new;
                tmass += (TMASS/noparticles);
              }
            }
          }
        }
      }
    }
  }
  cout << "tmass = " << tmass << endl;

  for(int i=0; i<20; i++){
    //cout << r[i] << endl;
  }
}


//fill opacity array
void opacity_fill(int NT, int NP, int NR){
  for(i=is(NT)+1; i<=ie(NT); i++){
    for(j=js(NP)+1; j<=je(NP); j++){
      for(k=ks(NR)+1; k<=ke(NR); k++){
        kappa[i][j][k] = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
      }
    }
  }
  //cout << "done opacity fill" << endl;
}

//calculate optical depth and fill array
void calculate_optical_depth_ana(int NT, int NP, int NR, double*** kappa, double*** d){
  for(i=is(NT)+1; i<=ie(NT); i++){
    for(j=js(NP)+1; j<=je(NP); j++){
      for(k=ks(NR)+1; k<=ke(NR); k++){
        t_ana[i][j][k]=t_ana[i][j][k-1]+(kappa[i][j][k]*d[i][j][k]*dRa[k]);
      }
    }
  }
}

void calculate_optical_depth_num(int NT, int NP, int NR, double*** kappa, double*** den){
  for(i=is(NT)+1; i<=ie(NT); i++){
    for(j=js(NP)+1; j<=je(NP); j++){
      for(k=ks(NR)+1; k<=ke(NR); k++){
        t_num[i][j][k]=t_num[i][j][k-1]+(kappa[i][j][k]*den[i][j][k]*dRa[k]);
      }
    }
  }
}
