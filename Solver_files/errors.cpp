#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"


using namespace std;

double delta( double value1, double value2){
    //evaluate error
    double del;
    del = fabs(value1 - value2);

    if (del < tol) {
	    return -1.0;
    } else {
	    return del / tol;
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
    vector <double> deltas(6);
    vector <double> order4(6), order5(6);

    //4th order
    order4 = new_variables(h, V, false);
    //5th order
    order5 = new_variables(h, V, true);

    for (unsigned int i = 0; i< 6; i++){
        deltas[i] = delta(order4[i], order5[i]);
    }
    return deltas;

}
