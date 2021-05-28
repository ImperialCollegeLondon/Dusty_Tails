#include <vector>
#include "spline.h"
#include <ostream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>
#include <stdlib.h> 
#include <chrono>
using namespace std;
random_device rd;
default_random_engine generator(rd()); // rd() provides a random seed
uniform_real_distribution<double> distribution_r(0.0,5.0);
uniform_real_distribution<double> distribution_t(0.0,4.0);
uniform_real_distribution<double> distribution_p(0.0,3.0);

int cells;
int cells_x = 4;
int cells_y = 4;

double cubicInterpolate ( vector <double> p, vector <double> s1, double x) {
    tk::spline s(s1,p, tk::spline::cspline);
    return s(x);
    
}

double bicubicInterpolate (vector< vector <double> > p, vector <double> s1, vector <double> s2, double x, double y) {
    vector <double> arr(s1.size());
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = cubicInterpolate(p[i], s2, y);
     }
    return cubicInterpolate(arr, s1, x);
}

double tricubicInterpolate (vector <vector< vector <double> > > p, vector <double> s1, vector <double> s2, vector <double> s3, double x, double y, double z) {
    vector <double> arr(s1.size());
     for (int i = 0; i < s1.size(); i++) {
         arr[i] = bicubicInterpolate(p[i], s2, s3, y, z);
     }
    return cubicInterpolate(arr, s1, x);
}


int main() {
    

    // default cubic spline (C^2) with natural boundary conditions (f''=0)
    vector< vector <vector <double>>> tau = {{{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16}, {17,18,19,20}}, {{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16}, {17,18,19,20}}, \
                                            {{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16}, {17,18,19,20}},{{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16},  {17,18,19,20}}, \
                                            {{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16},  {17,18,19,20}}, {{1,2,3,4}, {5,6,7,8}, {9,10,11, 12}, {13,14,15,16},  {17,18,19,20}} };
    vector <double> xs = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    vector <double> ys = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector <double> zs = {0.0, 1.0, 2.0, 3.0};
    
    //double test = bicubicInterpolate(tau2, xs, ys, 1.0, 0.0);
    //double test = tricubicInterpolate(tau, xs, ys, zs, 0.0, 4.0, 2.0);
    //cout << tau[2][2][1] << endl;
    //cout << test << endl;
    vector < vector <double> > points;

    for (int i = 0; i< 15000; i++) {
        double r = distribution_r(generator); 
        double t = distribution_t(generator);
        double p = distribution_p(generator);
        points.push_back({r,t,p});
    }
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i< 15000; i++) {
        vector <double> point = points[i];
        tricubicInterpolate(tau, xs, ys, zs, point[0], point[1], point[2]);
    }
    
    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;

    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}
