#include "CAV.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>

using namespace std;

vector <double> g1,g2;
double largest_element;

//cell width
vector <double> DR(double NR){
  double mu=1.;
  double sd=mu/3.;
  for(i=0.; i<=NR; i++){
    a=i/10.;
    gauss_new = (1./(sd*sqrt(2*M_PI))*exp(-pow(a-mu,2.)/(2.*pow(sd,2.))));
    g1.push_back(gauss_new);
    //cout << gauss_new << endl;
  }

  double largest_element = g1[0];
  for(int i = 1; i < g1.size(); i++){
    if(g1[i] > largest_element){
         largest_element = g1[i];
    }
  }

  for(i=0.; i<=NR; i++){
    a=i/10.;
    gauss_new = 0.1*((largest_element+0.1)-g1[i]);
    g2.push_back(gauss_new);
    cout << gauss_new << endl;
  }

  return g2;
}

int main(){
  DR(20.);
}
