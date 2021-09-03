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
uniform_real_distribution<double> distribution_r(0.0,100.0);
uniform_real_distribution<double> distribution_t(0.0,25.0);
uniform_real_distribution<double> distribution_p(0.0,150.0);

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

vector < vector < tk:: spline > > splines_phi (vector <vector< vector <double> > > taus, vector <double> radii, vector <double> thetas, vector <double> phis) {
    vector < vector < tk::spline >> phi_splines;
    tk::spline s; 
    for (int i = 0; i<radii.size(); i++){
        phi_splines.push_back({});
        for ( int j=0; j<thetas.size(); j++) {
            tk::spline s(phis, taus[i][j], tk::spline::cspline);
            phi_splines[i].push_back(s);
        }
    }

    return phi_splines;
}

vector <vector <double> > phi_spline_result(vector < vector < tk:: spline >> splines , vector<double> radii, vector <double> thetas, vector<double> phis, double phi) {
    vector <vector <double> > phi_spline_values;
    for (int i = 0; i<radii.size(); i++){
        phi_spline_values.push_back({});
        for ( int j=0; j<thetas.size(); j++) {
            phi_spline_values[i].push_back(splines[i][j](phi));
        }
    } 
    return phi_spline_values;
}
vector < tk:: spline > splines_theta ( vector <vector <double>> phi_splines_p, vector<double> radii, vector<double> thetas) {
     
    vector < tk::spline > r_splines;
    for (int i =0; i < radii.size(); i++){
        tk::spline s(thetas, phi_splines_p[i], tk::spline::cspline);
        r_splines.push_back(s);
    }

    return r_splines;
}

vector <double> theta_spline_result( vector <tk:: spline> splines, vector<double> radii, vector<double> thetas, vector<double> phis , double theta) {
    vector <double>  theta_spline_values;
    for (int i = 0; i<radii.size(); i++){
            theta_spline_values.push_back(splines[i](theta));
       
        }
    return theta_spline_values;
    } 
    

double tau_p (vector<double> theta_splines_p, vector<double> radii, vector<double> thetas, vector<double> phis, double radius){
    tk:: spline s(radii, theta_splines_p);
    
    return s(radius);
}


int main() {
    

    // default cubic spline (C^2) with natural boundary conditions (f''=0)
    vector< vector <vector <double>>> tau;
    vector <double> xs = {};
    vector <double> ys = {};
    vector <double> zs = {};

    double counter = 0.0;    
    for (int i = 0; i < 100; i++) {
        xs.push_back(counter);
        counter = counter + 1.0;
    }

    counter = 0.0;
    for (int i = 0; i < 25; i++) {
        ys.push_back(counter);
        counter = counter + 1.0;
    }
    counter = 0.0;

    for (int i = 0; i < 150; i++) {
        zs.push_back(counter);
        counter = counter + 1.0;
    }

    counter = 0.0;
    for (int i=0; i < 100; i++) {
        tau.push_back({});
        for (int j=0; j< 25; j++) {
            tau[i].push_back({});
            for (int k=0; k<150; k++) {
                tau[i][j].push_back(k);
            }
        }
    }
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

    //auto start = std::chrono::high_resolution_clock::now();
    
    vector< vector <tk::spline >> s_phi;

    s_phi = splines_phi (tau, xs, ys, zs);
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i< 15000; i++) {
        vector <double> point = points[i];
    
        vector < vector <double> > s_phi_values;
        s_phi_values = phi_spline_result(s_phi , xs, ys, zs, point[2]);

        vector < tk:: spline > s_theta;

        s_theta =  splines_theta ( s_phi_values, xs, ys);

        vector <double> s_theta_values;

        s_theta_values = theta_spline_result( s_theta, xs, ys, zs , point[1]);

        double result;
        result =  tau_p (s_theta_values, xs, ys, zs, point[0]);

    }
    
    
    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    
}



