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

//grid cell widths
vector <double> g,DR,x,inv;
double suminv,inv_new,sumDR,nr,C;
double mu;
double sd;

//cell width
vector <double> make_gauss(double nr){
  for(double i=0.; i<=(C*sd); i+=(1./nr)){
    x.push_back(i);
    gauss_new = exp(-pow(i-mu,2.)/(2.*pow(sd,2.)));
    g.push_back(gauss_new);
    //cout << i << " " << gauss_new << endl; / prints the gaussian vector elements
  }
  //cout << "done make gauss" << endl;
  return g;
  return x;
}

vector <double> find_inv(double B,vector<double> g){
  for(double i=0.; i<g.size(); i++){
    //cout << g[i] << endl; / prints the gaussian vector elements
    inv_new = 1.-(B*g[i]);
    if (inv_new == 0.){
      cout << "0 in DR vector" << endl; //so that DR isn't 0
    }
    inv.push_back(inv_new);
    //cout << inv[i] << endl;
  }
  //cout << "done find inv" << endl;
  suminv = accumulate(inv.begin(),inv.end(),0.); //calculates sum of inverse gaussian values
  //cout << suminv << endl; //prints the sum of the inverse gaussian values
  return inv;
}

vector <double> find_DR(double A, vector<double> inv){
  for(double i=0.; i<inv.size(); i++){
    //cout << g[i] << endl; // prints the gaussian vector elements
    gauss_new = A*inv[i];
    DR.push_back(gauss_new);
    //cout << gauss_new << endl; // prints the inverse gaussian vector elements
  }
  //cout << "done find DR" << endl;
  sumDR = accumulate(DR.begin(),DR.end(),0.);
  //cout << sumDR << endl; //prints the sum of the DR values
  return DR;
}

//starting indices
int is(double NR){
  return 0;
}

//ending indices
int ie(double NR){
  return int(NR);
}

//building the grid
void build_grid(double NR, vector<double> DR)
{
  for (i=is(NR)+1; i<=ie(NR); i++){ // defines the cell edges
    //cout << DR[i] << endl;
    Ra_new = Ra[i-1]+DR[i-1];
    Ra.push_back(Ra_new);
    //cout << Ra[i] << endl;
  }

  for (i=is(NR); i<=ie(NR)-1; i++){ // defines the cell centers
    Rb_new = Ra[i]+DR[i]/2.;
    Rb.push_back(Rb_new);
    //cout << Rb[i] << endl;
  }

  for (i=is(NR);i<=ie(NR)-1;i++){ // defines the width of a cell
    dRa_new = Ra[i+1]-Ra[i];
    dRa.push_back(dRa_new);
    //cout << i << " " << dRa[i] << endl;
  }

  for (i=is(NR)+1;i<=ie(NR)-1;i++){ // defines the width between two cell centers
    dRb_new = Rb[i]-Rb[i-1];
    dRb.push_back(dRb_new);
    //cout << dRb[i] << endl;
  }
  //cout << "done build grid" << endl;
}

//fill density array
void density_fill(double NR)
{
  mean = 1.0; //mean is(NR) at a, which is(NR) 1 in dimensionless units
  sd = 0.2*mean; //arbitrary

  // Constructing density array
  for(i=is(NR); i<=ie(NR)-1; i++){
    gauss_new = (density*pow(a,3.)/Mstar_kg)*exp(-pow((Rb[i])-mean,2.)/(2.*pow(sd,2.)));
    gauss.push_back(gauss_new);
    d.push_back(gauss_new);
    //cout << gauss_new << endl;
  }
  //cout << "done density fill" << endl;
}

//fill opacity array
void opacity_fill(double NR)
{
  //constructing opacity array
  for(i=is(NR); i<=ie(NR)-1; i++){
    k_new = (3./4.)*(1./density_bulk)*(1./1.e-6)*(Mstar_kg/pow(a,2.));
    k.push_back(k_new);
  }
  //cout << "done opacity fill" << endl;
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
