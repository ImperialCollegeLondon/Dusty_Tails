#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"

using namespace std;
//g++ RK_init.cpp solver.cpp microphysics.cpp maths.cpp errors.cpp kvalues.cpp forces.cpp to compile


//some declaration of variables

vector <double> V0; //initial variables vector

int main() {

    double init_vel;

    init_vel = pow((G_dim*(m_planet)/(0.1*r_h)), 0.5); //inertial frame

    //Define initial position in dimensionless units

    double x0 = planet_x + 0.1*r_h;
    double y0 = 0.0;
    double z0 = 0.0;

    //Define initial velocity in dimensionless units

    double xdot0 = 0.0;
    double ydot0 = 0.0;
    double zdot0 = 0.0;

    //initial variables vector

    V0.push_back(x0);
    V0.push_back(y0);
    V0.push_back(z0);
    V0.push_back(xdot0);
    V0.push_back(ydot0);
    V0.push_back(zdot0);

    RK_solver(h0, V0, 0.0);


}
