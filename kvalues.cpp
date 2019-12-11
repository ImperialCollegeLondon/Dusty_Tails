#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
//#include "RK_variables.h"

#include "functions.h"

using namespace std;

double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x;
double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y;
double k1_z, k2_z, k3_z, k4_z, k5_z, k6_z, k7_z;

double k1_xdot, k2_xdot, k3_xdot, k4_xdot, k5_xdot, k6_xdot, k7_xdot;
double k1_ydot, k2_ydot, k3_ydot, k4_ydot, k5_ydot, k6_ydot, k7_ydot;
double k1_zdot, k2_zdot, k3_zdot, k4_zdot, k5_zdot, k6_zdot, k7_zdot;


void k_values(double h, vector <double> V, bool order5, vector <double> star, \
              vector <double> planet){
    //function to obtain several k values of RK-DP method

    centrifugal(V[0], V[1], V[2], V[3], V[4], V[5]);
    coriolis(V[0], V[1], V[2], V[3], V[4], V[5]);
		rad_pressure(V[0], V[1], V[2], V[3], V[4], V[5]);
		pr_drag(V[0], V[1], V[2], V[3], V[4], V[5]);

    //k1 values
    k1_xdot = h* acceleration( V[0] - star_pos[0], V[0] - planet_pos[0],  \
			                         V[0],V[1], V[2], \
															 coriol_vector[0], centri_vector[0], \
															 rad_vector[0], pr_vector[0],
                               star, planet);

		k1_x = h* V[3];

    k1_ydot = h* acceleration( V[1] - star_pos[1], V[1] - planet_pos[1],   \
			                         V[0], V[1], V[2], \
															 coriol_vector[1], centri_vector[1], \
															 rad_vector[1], pr_vector[1], star, planet);
    k1_y = h* V[4];

    k1_zdot = h* acceleration( V[2] - star_pos[2], V[2] - planet_pos[2], \
			                         V[0], V[1], V[2], \
			                         coriol_vector[2], centri_vector[2], \
															 rad_vector[2], pr_vector[2], star, planet);

    k1_z = h* V[5];

    //k2 values

		centrifugal(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);

		coriolis(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);

		rad_pressure(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);

		pr_drag(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);


    k2_xdot = h*acceleration( V[0] - star_pos[0] + a21*k1_x,  \
			                        V[0] - planet_pos[0] + a21*k1_x, \
			                        V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			                        coriol_vector[0], centri_vector[0], \
														  rad_vector[0], pr_vector[0], star, planet);

    k2_x = h* (V[3] + a21*k1_xdot);

    k2_ydot = h*acceleration( V[1] - star_pos[1] + a21*k1_y, \
                              V[1] - planet_pos[1] + a21*k1_y, \
															V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
															coriol_vector[1], centri_vector[1], \
														  rad_vector[1], pr_vector[1], star, planet);

    k2_y = h* (V[4] + a21*k1_ydot);

    k2_zdot = h*acceleration( V[2] -star_pos[2] + a21*k1_z, \
			                        V[2] - planet_pos[2] + a21*k1_z, \
			                        V[0] + a21*k1_x, V[1] + a21*k1_y,V[2] + a21*k1_z, \
															coriol_vector[2], centri_vector[2], \
														  rad_vector[2], pr_vector[2], star, planet);

    k2_z = h* (V[5] + a21*k1_zdot);


    //k3 values

		centrifugal(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot );

		coriolis(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot);

		rad_pressure(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot);

		pr_drag(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot);

    k3_xdot = h* acceleration( V[0] - star_pos[0] + a31*k1_x + a32*k2_x, \
			                         V[0] - planet_pos[0] + a31*k1_x + a32*k2_x, \
                               V[0] + a31*k1_x + a32*k2_x, \
                               V[1] + a31*k1_y + a32*k2_y, \
                               V[2] + a31*k1_z + a32*k2_z, \
														   coriol_vector[0], centri_vector[0], \
                               rad_vector[0], pr_vector[0], star, planet);

    k3_x = h*(V[3] + a31*k1_xdot + a32*k2_xdot);

    k3_ydot = h* acceleration( V[1] - star_pos[1] + a31*k1_y + a32*k2_y, \
                               V[1] - planet_pos[1] + a31*k1_y + a32*k2_y, \
                               V[0] + a31*k1_x + a32*k2_x, \
                               V[1] + a31*k1_y + a32*k2_y, \
                               V[2] + a31*k1_z + a32*k2_z, \
														   coriol_vector[1], centri_vector[1], \
                               rad_vector[1], pr_vector[1], star, planet);

    k3_y = h*(V[4] + a31*k1_ydot + a32*k2_ydot);

    k3_zdot = h* acceleration( V[2] - star_pos[2] + a31*k1_z + a32*k2_z, \
                               V[2] - planet_pos[2] + a31*k1_z + a32*k2_z, \
                               V[0] + a31*k1_x + a32*k2_x, \
                               V[1] + a31*k1_y + a32*k2_y, \
                               V[2] + a31*k1_z + a32*k2_z, \
														   coriol_vector[2], centri_vector[2], \
                               rad_vector[2], pr_vector[2], star, planet);

    k3_z = h*(V[5] + a31*k1_zdot + a32*k2_zdot);



    //k4 values

		centrifugal(V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
		            V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
		            V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
							  V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
							  V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
							  V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

		coriolis(V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
		            V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
		            V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
							  V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
							  V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
							  V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

    rad_pressure(V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
		            V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
		            V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
							  V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
							  V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
							  V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

    pr_drag(V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
		            V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
		            V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
							  V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot, \
							  V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot, \
							  V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);




    k4_xdot = h* acceleration( V[0] - star_pos[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[0] - planet_pos[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[0], centri_vector[0], \
                               rad_vector[0], pr_vector[0], star, planet);

    k4_x = h* (V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot);

    k4_ydot = h* acceleration( V[1] - star_pos[1] + a41*k1_y + a42*k2_y + a43*k3_y, \
                               V[1] - planet_pos[1] + a41*k1_y + a42*k2_y + a43*k3_y, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[1], centri_vector[1], \
                               rad_vector[1], pr_vector[1], star, planet );

    k4_y = h* (V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot);

    k4_zdot = h* acceleration( V[2] - star_pos[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
                               V[2] - planet_pos[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[2], centri_vector[2], \
                               rad_vector[2], pr_vector[2], star, planet);

    k4_z = h* (V[5] + a41*k1_zdot + a42*k2_zdot + a43*k3_zdot);

    //k5 values

		centrifugal(V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
		            V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
		            V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
							  V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
							  V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
							  V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

		coriolis(V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
		         V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
		         V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
						 V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
						 V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
						 V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

    rad_pressure(V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
		         V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
		         V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
						 V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
						 V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
						 V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

    pr_drag(V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
		         V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
		         V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
						 V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot, \
						 V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot, \
						 V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);

    k5_xdot = h* acceleration( V[0] - star_pos[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[0] - planet_pos[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[0], centri_vector[0], \
                               rad_vector[0], pr_vector[0], star, planet);

    k5_x = h* (V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot);

    k5_ydot = h* acceleration( V[1] - star_pos[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y, \
                               V[1] - planet_pos[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[1], centri_vector[1], \
                               rad_vector[1], pr_vector[1], star, planet);

    k5_y = h* (V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot);

    k5_zdot = h* acceleration( V[2] - star_pos[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
                               V[2] - planet_pos[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[2], centri_vector[2], \
                               rad_vector[2], pr_vector[2], star, planet);

    k5_z = h* (V[5] + a51*k1_zdot + a52*k2_zdot + a53*k3_zdot + a54*k4_zdot);


    //k6 values

		centrifugal(V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
		            V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
		            V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
							  V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
							  V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
							  V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);

		coriolis(V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
		         V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
		         V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
						 V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
						 V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
						 V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);


    rad_pressure(V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
		         V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
		         V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
						 V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
						 V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
						 V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);

    pr_drag(V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
		         V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
		         V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
						 V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot, \
						 V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot, \
						 V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);



    k6_xdot = h* acceleration( V[0] - star_pos[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[0] - planet_pos[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[0], centri_vector[0], \
                               rad_vector[0], pr_vector[0], star, planet);

    k6_x = h* (V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot);

    k6_ydot = h* acceleration( V[1] - star_pos[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[1] - planet_pos[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[1], centri_vector[1], \
                               rad_vector[1], pr_vector[1], star, planet);

    k6_y = h* (V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot);

    k6_zdot = h* acceleration( V[2] - star_pos[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
                               V[2] - planet_pos[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[2], centri_vector[2], \
                               rad_vector[2], pr_vector[2], star, planet);

    k6_z = h* (V[5] + a61*k1_zdot + a62*k2_zdot + a63*k3_zdot + a64*k4_zdot + a65*k5_zdot);


    if (order5 == true) {

			centrifugal(V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
			            V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
			            V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
								  V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot, \
								  V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot, \
								  V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);

			coriolis(V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
			         V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
			         V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
							 V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot, \
							 V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot, \
							 V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);

      rad_pressure(V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
			         V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
			         V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
							 V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot, \
							 V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot, \
							 V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);

      pr_drag(V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
			         V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
			         V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
							 V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot, \
							 V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot, \
							 V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);





      k7_xdot = h* acceleration( V[0] - star_pos[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                 V[0] - planet_pos[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                 V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                 V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                 V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																 coriol_vector[0], centri_vector[0], \
                                 rad_vector[0], pr_vector[0], star, planet);

      k7_x = h* (V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot);

      k7_ydot = h* acceleration( V[1] - star_pos[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                 V[1] - planet_pos[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                 V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                 V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                 V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																 coriol_vector[1], centri_vector[1], \
                                 rad_vector[1], pr_vector[1], star, planet);

        k7_y = h* (V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot);

        k7_zdot = h* acceleration( V[2] - star_pos[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
                                   V[2] - planet_pos[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																   coriol_vector[2], centri_vector[2], \
                                   rad_vector[2], pr_vector[2], star, planet);


        k7_z = h* (V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);

    }
}
