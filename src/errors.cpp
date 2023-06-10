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
#include <omp.h>


//calculates errors on numerical solver to accept,increase 
//or decrease current "small" timestep
using namespace std;
double error( vector<double> v1, vector<double> v2, bool print) {
  double err, frac;
  double oom1, oom2, oom;
  //double tol;
  vector <double> tol(7);
  tol[0] = 1.0e-5;
  tol[1] = 1.0e-8;
  tol[2] = 1.0e-8;
  tol[3] = 1.0e-8;
  tol[4] = 1.0e-8;
  tol[5] = 1.0e-8;
  tol[6] = 1.0e-8;
  err = 0.0;
  for (unsigned int i = 0; i < 7; i++){
    // oom1 = log10(abs(v1[i]));
    // //cout << "oom1 " << oom1 << endl;
    // oom2 = log10(abs(v2[i]));
    // //cout << "oom2 " << oom2 << endl;
    // oom = min(oom1, oom2);
    // oom = round(oom);
    // //cout << "oom " << oom << endl;
    // tol = pow(10.0, oom) * 1.0e-4 ;
    // if (isnan(tol)){
    //   cout << "tol is nan " << endl;
    //   abort();
    // }
    // if (i==6) {
    //   tol = 1.0e-8;
    // }
    
    //cout << "tol " << tol << endl;
    frac = fabs(v1[i] - v2[i]) / (tol[i]+max(abs(v1[i]), abs(v2[i])) * tol[i]);
    err = err + pow( (1./7.)*pow(frac,2.0),0.5);
    // if (print)
    // { 
    //   cout << "variable " << i << endl;
    //   //cout << "tol " << tol[i] << endl;
    //   cout << "value 1 " << v1[i] << endl;
    //   cout << "value 2 " << v2[i] << endl;
    //   cout << "error " << frac << endl;
    // } 
  }

  
  //cout << "error " << err << endl;
  return err;
}


tuple<double, int> error_max(double h, vector <double> V){
    vector <double> errors(7);
    vector <double> order4(7), order5(7);

    
    double err;
    int max_index;
    bool debug, print;
    debug = false;
    //4th order
    order4 = new_variables(h, V, false, false);
    //5th order
    order5 = new_variables(h, V, true, false);
    if (h<1.0e-8) {
      print = true;
    } else {
      print = false;
    }
    err = error(order4, order5, print);
    max_index = 0;
    return make_tuple(err, max_index);

}

vector <double> new_step_size(tuple<double,int> errors, double h, 
                              bool reject, vector <double> V, double err_old, int counter){
  double h_new;
  double err;
  double alpha, beta, safe;
  double minscale, maxscale, scale;
  
  minscale = 0.5;
  maxscale = 10.0;
  safe = 0.90;
  beta = 0.08;
  alpha = 0.2 - beta*0.75;
  vector <double> steps;
  err = get<0>(errors);
  
  // if (counter>5) {
  //   cout << "oops struggling in errors " << endl;
  //   cout << "old error " << err_old << endl;
  //   cout << "current error " << err << endl;
  //   cout << "stepsize " << h << endl;
  // }

  if (isnan(err)){
    cout << "max err is NaN " << endl;
    abort();
  }

  if (err <= 1.0) {
    //sucessful step choice
    
    if (err == 0.0) {
        scale = maxscale;
        
    } 
    else{
        scale = safe * pow(err, -alpha) * pow(err_old, beta);
        if (scale < minscale) {scale = minscale;}
        if (scale > maxscale) {scale = maxscale;}
        if (isnan(scale))  {scale = maxscale;}
        //scale = maxscale;
    }
    if (reject) {
        h_new = h * min(scale, 1.0);
    }else{
        h_new = h * scale;
        err_old = max(err, 1.0e-4);
    }
    // if (h < 1.0e-8)
    // {
    //     cout << "sucess!" << endl;
    //     cout << "time step is quite small " << h << endl;
    //     cout << "new time step " << h_new << endl;
    //     cout << "scale " << scale << endl;
    //     //cout << "old error " << err_old << endl;
    //     cout << "new error " << err << endl;
    //     cout << "reject " << reject << endl;
    //     cout << "counter " << counter << endl;
    // }
    // if (counter > 5)
    // {
    //   cout << "scale " << scale << endl;
    //   cout << "new h " << h_new << endl;
    // }  
        if (isinf(h_new) or isnan(h_new))
        {
          cout << "new time step overflows!" << h_new << endl;
          cout << " scale " << scale << endl;
          cout << " old time step " << h << endl;
          cout << "error estimated " << err << endl;
          cout << "old error " << err_old << endl;
          cout << "aborting " << endl;
          abort();
        }
        
        steps = {h, h_new, err};
        return steps;
    }
    else
    {
        // failed
        //  cout << " failed, error is " << err << endl;
        //  cout << "old error was " << err_old << endl;
        //  cout << "old step was " << h << endl;
        //  cout <<  " scale " <<  safe * pow(err, -alpha) << endl;
        counter = counter + 1;
        scale = max(safe * pow(err, -alpha), minscale);
        // if (h < 1.0e-8)
        // {
        //   cout << "time step is quite small " << h << endl;
        //   cout << "old h " << h / scale << endl;
        //   cout << "scale " << scale << endl;
        //   cout << "old error " << err << endl;
        //   cout << "new error " << get<0>(errors) << endl;
        //   cout << "reject " << reject << endl;
        //   cout << "counter " << counter << endl;
          
        // }
        h = h * scale;
        
        // cout << "new h will be " << h << endl;
        errors = error_max(h, V);
        
        return new_step_size(errors, h, true, V, err, counter);
    }
}
