#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "particle.h"
#include "functions.h"

using namespace std;

vector <double> to_vector(double x, double y, double z){
  //function to put any 3 values into a vector
  vector <double> some_vector;
	some_vector.clear();

	some_vector.push_back(x);
	some_vector.push_back(y);
	some_vector.push_back(z);

	return some_vector;
}

double scalar(double x, double y, double z){
        double s;
        s = pow( pow(x, 2.) + pow(y, 2.) + pow(z, 2.), 0.5 );
        return s;
}

vector <double> cross_product(double m1, double m2, double m3, double n1, double n2, double n3){
  //function to calculate the cross product between two vectors
  vector <double> cp_vector;
  cp_vector.clear();

	double i_new, j_new, k_new;

	i_new = m2*n3 - m3*n2;
	j_new = m3*n1 - m1*n3;
	k_new = m1*n2 - m2*n1;

	cp_vector.push_back(i_new);
	cp_vector.push_back(j_new);
	cp_vector.push_back(k_new);

	return cp_vector;


}

double dot_product(vector <double> n,  vector <double> m){
     //function to calculate the dot product between two vectors
     double dp = (n[0]*m[0]) + (n[1]*m[1]) + (n[2]*m[2]);

		 return dp;

}
