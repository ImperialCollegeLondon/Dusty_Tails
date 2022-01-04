#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include "opacities.h"
#include <tuple>


//calculates errors on numerical solver to accept,increase 
//or decrease current "small" timestep
using namespace std;

double error( double value1, double value2){
    //evaluate relative error
    double err;
   
    err = fabs(value1 - value2)/abs(value1);
    
    return err;
}

tuple<double, int> error_max(double h, vector <double> V){
    vector <double> errors(7);
    vector <double> order4(7), order5(7);
    double err_max;
    int max_index;

    //4th order
    order4 = new_variables(h, V, false);
    //5th order
    order5 = new_variables(h, V, true);

    for (unsigned int i = 0; i < 7; i++){
        
        errors[i] = error(order4[i], order5[i]);

    }

    err_max = *max_element(errors.begin(), errors.end());
    max_index = max_element(errors.begin(),errors.end()) - errors.begin();
    
    return make_tuple(err_max, max_index);

}

vector <double> new_step_size(tuple<double,int> errors, double h_old, int fail_status, vector <double> V){
  double rho;
  double h_new;
  double max_err;
  double tol1;
 
  vector <double> steps;
  max_err = get<0>(errors);
  //cout << max_err << "    " << get<1>(errors) << endl;
  //if (get<1>(errors)<6) {
    //tol1 = 1.0e-4;
 // } else {
    //tol1 = 1.0e-7;
  //}
  tol1 = 1.0e-3;
  rho = 1.25 * pow((max_err / tol1), 1.0/5.0);

  if (isnan(max_err)){
    cout << "max err is NaN " << endl;
    abort();
  }

  if (max_err <= tol1){
    if (rho > 0.2){
      h_new = h_old / rho;
      steps = {h_old, h_new};
      return steps;

    } else {
  
      h_new = 5.0 * h_old;
      steps = {h_old, h_new};
      return steps;
    }
  } else {
    if (fail_status == 0){
      h_new = max(0.1, 1./rho) * h_old;
      errors = error_max(h_new, V);
      return new_step_size(errors, h_new, 1, V);
    } else {
      h_new = 0.5 * h_old;
      errors = error_max(h_new, V);
      return new_step_size(errors, h_new, 1, V);
    }
  }

}
