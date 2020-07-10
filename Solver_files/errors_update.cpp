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

double error( double value1, double value2){
    //evaluate error
    double err;
    err = fabs(value1 - value2);
    //cout << "error " << err/tol << endl;
    return err / tol;
}

double error_max(double h, vector <double> V){
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

    return err_max;

}

vector <double> new_step_size(double max_err, double h_old, int fail_status, vector <double> V){
  double rho;
  double h_new;
  double new_max_err;
  vector <double> steps;


  //cout << "fail " << fail_status << endl;

  rho = 1.25 * pow((max_err / tol), 1.0/5.0);



  if (isnan(max_err)){
    return {-1.0, -1.0};
  }

  if (max_err <= tol){
    if (rho > 0.2){
      h_new = h_old / rho;
      //h_new = 0.001;
      steps = {h_old, h_new};
      return steps;
    } else {
      //h_new = 0.001;
      h_new = 5.0 * h_old;
      steps = {h_old, h_new};
      return steps;
    }
  } else {
    if (fail_status == 0){
      h_new = max(0.1, 1./rho) * h_old;
      new_max_err = error_max(h_new, V);
      return new_step_size(new_max_err, h_new, 1, V);
    } else {
      h_new = 0.5 * h_old;
      new_max_err = error_max(h_new, V);
      return new_step_size(new_max_err, h_new, 1, V);
    }
  }

}
