#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"

using namespace std;

double mu = 101.961; //for corundum, molecular weight
double alpha = 0.1;
double A = 77365.;
double B = 39.3;

const double kb = 1.381e-16;
const double amu = 1.661e-24;

double rho = 4.0;

vector <double> arange( double start, double stop, double step) {
    vector <double> values;
    for (double value = start; value < stop; value += step){
        values.push_back(value);
    }
    return values;
}

double timescale(double T){
    double J;
    J = -(alpha/rho)* exp((-A/T) + B) * pow((mu * amu)/(2.0*PI*kb*T), 0.5);
    return size/J;
}

int main(){

    vector <double> temperatures;
    vector <double> times;
    temperatures = arange( 1200.0, 8000.0, 5.);

    ofstream file("sublimation_data.txt");
    vector <double>::iterator it;

    for (it = temperatures.begin(); it != temperatures.end(); it++){
        double t_now;
        t_now = timescale(*it)/T;
        times.push_back(fabs(t_now));
        file << *it << ",";
        file << fabs(t_now) << "\n";

    }

}
