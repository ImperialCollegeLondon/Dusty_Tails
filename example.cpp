#include "CAV.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>

using namespace std;

vector <double> g,DR,x;
double largest_element;
double mu = 1.;
double sd = Rmax/6.;

//cell width
vector <double> make_gauss(double NR){
  for(double i=0.; i<=(6.*sd); i+=(1./NR)){
    x.push_back(i);
    gauss_new = exp(-pow(i-mu,2.)/(2.*pow(sd,2.)));
    g.push_back(gauss_new);
    //cout << i << " " << gauss_new << endl; / prints the gaussian vector elements
  }
  return g;
  return x;
}

vector <double> find_DR(double A,double B, vector<double> g){
  for(double i=0.; i<g.size(); i++){
    //cout << g[i] << endl; / prints the gaussian vector elements
    gauss_new = A*(1.-g[i])+B;
    DR.push_back(gauss_new);
    //cout << gauss_new << endl; / prints the inverse gaussian vector elements
  }
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

void file_creator_inv(vector<double> DR, vector <double> x) {//textfile of density vs Rb
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


int main(){
  make_gauss(NR=100.);
  find_DR(A=10.,B=0.1,g);
  file_creator_gaussian(g,x);
  file_creator_inv(DR,x);

  return 0;
}
