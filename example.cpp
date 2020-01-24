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

vector <double> g,DR,x,inv;
double mu = 1.;
double sd = Rmax/20.;
double suminv,inv_new,sumDR;
double nr;

//cell width
vector <double> make_gauss(double nr){
  for(double i=0.; i<=(20.*sd); i+=(1./nr)){
    x.push_back(i);
    gauss_new = exp(-pow(i-mu,2.)/(2.*pow(sd,2.)));
    g.push_back(gauss_new);
    //cout << i << " " << gauss_new << endl; / prints the gaussian vector elements
  }
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
  suminv = accumulate(inv.begin(),inv.end(),0.); //calculates sum of inverse gaussian values
  //cout << suminv << endl; //prints the sum of the inverse gaussian values
  return inv;
}

vector <double> find_DR(double A, vector<double> inv){
  for(double i=0.; i<inv.size(); i++){
    //cout << g[i] << endl; // prints the gaussian vector elements
    //cout << A << endl;
    gauss_new = A*inv[i];
    DR.push_back(gauss_new);
    //cout << gauss_new << endl; // prints the inverse gaussian vector elements
  }
  sumDR = accumulate(DR.begin(),DR.end(),0.);
  //cout << sumDR << endl; //prints the sum of the DR values
  return DR;
}

void file_creator_gaussian(vector<double> g, vector <double> x) {//textfile of density vs Rb
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

void file_creator_DR(vector<double> DR, vector <double> x) {//textfile of density vs Rb
  ofstream myfile ("DR.txt");
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


int main(){
  make_gauss(nr=20.);
  find_inv(B=0.9,g);
  find_DR(A=(2./suminv),inv);
  file_creator_gaussian(g,x);
  file_creator_DR(DR,x);

  return 0;
}
