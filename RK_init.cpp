#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
//#include "RK_variables.h""
#include "functions.h"

using namespace std;

//some declaration of variables


vector <double> V0; //initial variables vector
vector <double> star_pos, planet_pos; //star and planet position in CoM frame
vector <double> cp_vector, centri_vector, coriol_vector, vrad, rad_vector, s_unit, pr_vector;
vector <double> x_positions, y_positions, z_positions;

int main() {

    double init_vel;

    init_vel = pow((G_dim*(m_planet)/(0.1*r_h)), 0.5); //inertial frame

    //Define initial position in dimensionless units

		double star_x  = -((m_planet) / (m_planet + 1.0));
		double star_y = 0.0;
		double star_z = 0.0;

		double planet_x = 1.0 / (m_planet + 1.0);
		double planet_y = 0.0;
		double planet_z = 0.0;

		star_pos.push_back(star_x);
		star_pos.push_back(star_y);
		star_pos.push_back(star_z);

		planet_pos.push_back(planet_x);
		planet_pos.push_back(planet_y);
		planet_pos.push_back(planet_z);

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

    RK_solver(h0, V0, 0.0, star_pos, planet_pos);


}
