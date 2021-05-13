#include <vector>
#include "spline.h"
#include <ostream>
#include <iostream>
using namespace std;

int main() {
    vector<double> X, Y;

    // default cubic spline (C^2) with natural boundary conditions (f''=0)

    X = {1., 2., 3., 4., 5.};
    Y = {2., 4., 6., 8., 10.};

    tk::spline s;			// X needs to be strictly increasing
        
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::first_deriv, 1.0);
    s.set_points(X,Y);
    s.make_monotonic();
    double value=s(1.3);		// interpolated value at 1.3
    double deriv=s.deriv(1,1.3);	// 1st order derivative at 1.3
    vector<double> solutions = s.solve(3.0);	// solves s(x)=0.0

    cout << solutions[0] << endl;

}