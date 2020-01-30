#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
//#include "RK_variables.h"
#include "functions.h"
#include "particle.h"


using namespace std;

double delta( double value1, double value2){
    //evaluate error
    double del;
    del = fabs(value1 - value2);

    if (del < (tol + max(fabs(value1), fabs(value2))*tol)) {
	    return -1.0;
    } else {
	    return del / (tol + max(fabs(value1), fabs(value2))*tol);
    }

}

double h_optimal(vector <double> deltas_list, double h){
    //function to evaluate optimal h
    double max_delta, h_opt;

    max_delta = *max_element(deltas_list.begin(), deltas_list.end()); //includes tolerance

    if (max_delta == -1.0){
	    return -1.0;
    } else {

	    h_opt = S * h * pow(1/max_delta, 1.0/5.0);
	    return h_opt;
    }


}

vector <double> h_check(double h, vector <double> V){
    vector <double> deltas;
    vector <double> order4, order5;

    double x4, y4, z4, xdot4, ydot4, zdot4;
    double x5, y5, z5, xdot5, ydot5, zdot5;

    double xdot_err, ydot_err, zdot_err;
    double x_err, y_err, z_err;

    deltas.clear();

    //4th order
    order4 = new_variables(h, V, false);

    x4 = order4[0];
    y4 = order4[1];
    z4 = order4[2];

    xdot4 = order4[3];
    ydot4 = order4[4];
    zdot4 = order4[5];

    //5th order
    order5 = new_variables(h, V, true);

    x5 = order5[0];
    y5 = order5[1];
    z5 = order5[2];

    xdot5 = order5[3];
    ydot5 = order5[4];
    zdot5 = order5[5];

    //error on xdot
    xdot_err = delta(xdot4, xdot5);
    deltas.push_back(xdot_err);

    //error on ydot
    ydot_err = delta(ydot4, ydot5);
    deltas.push_back(ydot_err);

    //error on zdot
    zdot_err = delta(zdot4, zdot5);
    deltas.push_back(zdot_err);


    //error on x
    x_err = delta(x4, x5);
    deltas.push_back(x_err);

    //error on y
    y_err = delta(y4, y5);
    deltas.push_back(y_err);

    //error on z
    z_err = delta(z4, z5);
    deltas.push_back(z_err);

    return deltas;

}
