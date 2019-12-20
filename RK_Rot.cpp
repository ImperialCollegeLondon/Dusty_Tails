#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "RK_variables.h"

using namespace std;

//some declaration of variables

vector <double> V0; //initial variables vector
vector <double> dp_vector, centri_vector, coriol_vector;
vector <double> x_positions, y_positions, z_positions;

vector <double> x, y, z, V_new;

vector <double> deltas, delta_values;
vector <double> time_plot, semis, timing, as;

double scalar(double x, double y, double z);

vector <double> cross_product(double m1, double m2, double m3, double n1, double n2, double n3);

//end of variables declaration

vector <double> cross_product(double m1, double m2, double m3, double n1, double n2, double n3){

	dp_vector.clear();

	double i_new, j_new, k_new;

	i_new = m2*n3 - m3*n2;
	j_new = m3*n1 - m1*n3;
	k_new = m1*n2 - m2*n1;

	dp_vector.push_back(i_new);
	dp_vector.push_back(j_new);
	dp_vector.push_back(k_new);

	return dp_vector;


}



double semimajor(double period) { //function to evaluate semi-major axis value, period in days, a in AU

   semi = pow(7.496e-6 * Mstar_sun * pow(Period_days, 2.0), (1.0/3.0)); //in AU
   return semi*AU_to_m; //output in meters
}


vector <double> centrifugal(double x, double y, double z, double vx, double vy, double vz){

  vector <double> omxr, finalx;

	centri_vector.clear();

	omxr = cross_product(0.0, 0.0, 2.0*PI, x, y, z);

  finalx = cross_product(0.0, 0.0, 2.0*PI, omxr[0], omxr[1], omxr[2]);

	centri_vector.push_back(finalx[0]);
	centri_vector.push_back(finalx[1]);
	centri_vector.push_back(finalx[2]);

	return centri_vector;
}


vector <double> coriolis(double x, double y, double z, double vx, double vy, double vz){
  vector <double> omxv,  rvxr, finalx;

	coriol_vector.clear();
  omxv = cross_product(0.0, 0.0, 2.0*PI, vx, vy, vz);

	coriol_vector.push_back(omxv[0]);
	coriol_vector.push_back(omxv[1]);
	coriol_vector.push_back(omxv[2]);

	return coriol_vector;

}


double acceleration( double pos, double x, double y, double z, double coriol, double centri){

  vel_dot = ((-G_dim * pos * (1 - beta)) / pow( ( pow(x, 2.0) + pow(y, 2.0) + \
	            pow(z, 2.0) ) , (3.0/2.0) )) - centri - 2.0*coriol;
  return vel_dot;
}


void k_values(double h, vector <double> V, bool order5){


    centrifugal(V[0], V[1], V[2], V[3], V[4], V[5]);
    coriolis(V[0], V[1], V[2], V[3], V[4], V[5]);

    //k1 values
    k1_xdot = h* acceleration( V[0], V[0], V[1], V[2], coriol_vector[0], \
			centri_vector[0]);

		k1_x = h* V[3];

    k1_ydot = h* acceleration( V[1], V[0], V[1], V[2], coriol_vector[1], \
			centri_vector[1]);
    k1_y = h* V[4];

    k1_zdot = h* acceleration( V[2], V[0], V[1], V[2], coriol_vector[2], \
			centri_vector[2]);
    k1_z = h* V[5];

    //k2 values

		centrifugal(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);

		coriolis(V[0] + a21*k1_x, V[1] + a21*k1_y, V[2] + a21*k1_z, \
			V[3] + a21*k1_xdot, V[4] + a21*k1_ydot, V[5] + a21*k1_zdot);

    k2_xdot = h*acceleration( V[0] + a21*k1_x, V[0] + a21*k1_x, V[1] + a21*k1_y,\
			 V[2] + a21*k1_z, coriol_vector[0], centri_vector[0]);

    k2_x = h* (V[3] + a21*k1_xdot);

    k2_ydot = h*acceleration( V[1] + a21*k1_y, V[0] + a21*k1_x, V[1] + a21*k1_y, \
			V[2] + a21*k1_z, coriol_vector[1], centri_vector[1]);

    k2_y = h* (V[4] + a21*k1_ydot);

    k2_zdot = h*acceleration( V[2] + a21*k1_z, V[0] + a21*k1_x, V[1] + a21*k1_y,\
			 V[2] + a21*k1_z, coriol_vector[2], centri_vector[2]);

    k2_z = h* (V[5] + a21*k1_zdot);


    //k3 values

		centrifugal(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot );

		coriolis(V[0] + a31*k1_x + a32*k2_x, V[1] + a31*k1_y + a32*k2_y, \
			V[2] + a31*k1_z + a32*k2_z, V[3] + a31*k1_xdot + a32*k2_xdot, \
			V[4] + a31*k1_ydot + a32*k2_ydot, V[5] + a31*k1_zdot + a32*k2_zdot);

    k3_xdot = h* acceleration( V[0] + a31*k1_x + a32*k2_x, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z, \
														coriol_vector[0], centri_vector[0]);

    k3_x = h*(V[3] + a31*k1_xdot + a32*k2_xdot);

    k3_ydot = h* acceleration( V[1] + a31*k1_y + a32*k2_y, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z, \
														  coriol_vector[1], centri_vector[1]);

    k3_y = h*(V[4] + a31*k1_ydot + a32*k2_ydot);

    k3_zdot = h* acceleration( V[2] + a31*k1_z + a32*k2_z, V[0] + a31*k1_x + a32*k2_x, \
                              V[1] + a31*k1_y + a32*k2_y, V[2] + a31*k1_z + a32*k2_z, \
														  coriol_vector[2], centri_vector[2]);

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

    k4_xdot = h* acceleration( V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[0], centri_vector[0]);

    k4_x = h* (V[3] + a41*k1_xdot + a42*k2_xdot + a43*k3_xdot);

    k4_ydot = h* acceleration( V[1] + a41*k1_y + a42*k2_y + a43*k3_y, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[1], centri_vector[1]);

    k4_y = h* (V[4] + a41*k1_ydot + a42*k2_ydot + a43*k3_ydot);

    k4_zdot = h* acceleration( V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
                               V[0] + a41*k1_x + a42*k2_x + a43*k3_x, \
                               V[1] + a41*k1_y + a41*k2_y + a43*k3_y, \
                               V[2] + a41*k1_z + a42*k2_z + a43*k3_z, \
														   coriol_vector[2], centri_vector[2]);

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

    k5_xdot = h* acceleration( V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[0], centri_vector[0]);

    k5_x = h* (V[3] + a51*k1_xdot + a52*k2_xdot + a53*k3_xdot + a54*k4_xdot);

    k5_ydot = h* acceleration( V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[1], centri_vector[1]);

    k5_y = h* (V[4] + a51*k1_ydot + a52*k2_ydot + a53*k3_ydot + a54*k4_ydot);

    k5_zdot = h* acceleration( V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
                               V[0] + a51*k1_x + a52*k2_x + a53*k3_x + a54* k4_x, \
                               V[1] + a51*k1_y + a52*k2_y + a53*k3_y + a54* k4_y,\
                               V[2] + a51*k1_z + a52*k2_z + a53*k3_z + a54* k4_z, \
														   coriol_vector[2], centri_vector[2]);

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

    k6_xdot = h* acceleration( V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[0], centri_vector[0]);

    k6_x = h* (V[3] + a61*k1_xdot + a62*k2_xdot + a63*k3_xdot + a64*k4_xdot + a65*k5_xdot);

    k6_ydot = h* acceleration( V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[1], centri_vector[1]);

    k6_y = h* (V[4] + a61*k1_ydot + a62*k2_ydot + a63*k3_ydot + a64*k4_ydot + a65*k5_ydot);

    k6_zdot = h* acceleration( V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
                               V[0] + a61*k1_x + a62*k2_x + a63*k3_x + a64*k4_x + a65*k5_x, \
                               V[1] + a61*k1_y + a62*k2_y + a63*k3_y + a64*k4_y + a65*k5_y, \
                               V[2] + a61*k1_z + a62*k2_z + a63*k3_z + a64*k4_z + a65*k5_z, \
														   coriol_vector[2], centri_vector[2]);

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

        k7_xdot = h* acceleration( V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																   coriol_vector[0], centri_vector[0]);

        k7_x = h* (V[3] + a71*k1_xdot + a73*k3_xdot + a74*k4_xdot + a75* k5_xdot + a76*k6_xdot);

        k7_ydot = h* acceleration( V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																   coriol_vector[1], centri_vector[1]);

        k7_y = h* (V[4] + a71*k1_ydot + a73*k3_ydot + a74*k4_ydot + a75* k5_ydot + a76*k6_ydot);

        k7_zdot = h* acceleration( V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
                                   V[0] + a71*k1_x + a73*k3_x + a74*k4_x + a75*k5_x + a76*k6_x, \
                                   V[1] + a71*k1_y + a73*k3_y + a74*k4_y + a75*k5_y + a76*k6_y, \
                                   V[2] + a71*k1_z + a73*k3_z + a74*k4_z + a75*k5_z + a76*k6_z, \
																   coriol_vector[2], centri_vector[2]);


        k7_z = h* (V[5] + a71*k1_zdot + a73*k3_zdot + a74*k4_zdot + a75* k5_zdot + a76*k6_zdot);

    }
}

double new_variables(double h, vector <double> V, bool order5){

    if (order5 == false) {

     k_values(h, V, false);

     x_new = V[0] + b1*k1_x + b3* k3_x + b4* k4_x + b5* k5_x + b6* k6_x;
     y_new = V[1] + b1*k1_y + b3* k3_y + b4* k4_y + b5* k5_y + b6* k6_y;
     z_new = V[2] + b1*k1_z + b3* k3_z + b4* k4_z + b5* k5_z + b6* k6_z;

     xdot_new = V[3] + b1*k1_xdot  + b3* k3_xdot + b4* k4_xdot + b5* k5_xdot + b6* k6_xdot;
     ydot_new = V[4] + b1*k1_ydot  + b3* k3_ydot + b4* k4_ydot + b5* k5_ydot + b6* k6_ydot;
     zdot_new = V[5] + b1*k1_zdot  + b3* k3_zdot + b4* k4_zdot + b5* k5_zdot + b6* k6_zdot;

    } else {

        k_values(h, V, true);

        x_new = V[0] + bs1*k1_x + bs3*k3_x + bs4*k4_x + bs5*k5_x + bs6*k6_x + bs7*k7_x;
        y_new = V[1] + bs1*k1_y + bs3*k3_y + bs4*k4_y + bs5*k5_y + bs6*k6_y + bs7*k7_y;
        z_new = V[2] + bs1*k1_z + bs3*k3_z + bs4*k4_z + bs5*k5_z + bs6*k6_z + bs7*k7_z;

        xdot_new = V[3] + bs1*k1_xdot  + bs3* k3_xdot + bs4* k4_xdot + bs5* k5_xdot + bs6* k6_xdot + bs7*k7_xdot;
        ydot_new = V[4] + bs1*k1_ydot  + bs3* k3_ydot + bs4* k4_ydot + bs5* k5_ydot + bs6* k6_ydot + bs7*k7_ydot;
        zdot_new = V[5] + bs1*k1_zdot  + bs3* k3_zdot + bs4* k4_zdot + bs5* k5_zdot + bs6* k6_zdot + bs7*k7_zdot;


    }

    return x_new, y_new, z_new, xdot_new, ydot_new, zdot_new;
}



double delta( double value1, double value2){

    del = fabs(value1 - value2);

    if (del < (tol + fabs(value1)*tol) ) {
	    return -1.0;
    } else {
	return del / (tol + fabs(value1)*tol);
    }

}

double h_optimal(vector <double> deltas_list, double h){


    max_delta = *max_element(deltas_list.begin(), deltas_list.end()); //includes tolerance

    if (max_delta == -1.0){
	    return h;
    } else {

	    h_opt = S * h * pow(1/max_delta, 1.0/5.0);
	    return h_opt;
    }


}

vector <double> h_check(double h, vector <double> V){

    deltas.clear();

    //4th order
    new_variables(h, V, false);

    x4 = x_new;
    y4 = y_new;
    z4 = z_new;

    xdot4 = xdot_new;
    ydot4 = ydot_new;
    zdot4 = zdot_new;

    //5th order
    new_variables(h, V, true);

    x5 = x_new;
    y5 = y_new;
    z5 = z_new;
    xdot5 = xdot_new;
    ydot5 = ydot_new;
    zdot5 = zdot_new;


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


vector <double> next_step(double h, vector <double> V) {

    //ensure vector is  clear
    V_new.clear();

    new_variables(h, V, false);

    V_new.push_back(x_new);
    V_new.push_back(y_new);
    V_new.push_back(z_new);

    V_new.push_back(xdot_new);
    V_new.push_back(ydot_new);
    V_new.push_back(zdot_new);

    return V_new;
}


void RK_solver(double h0, vector <double> V_0, double t_0){

    double t = t_0;
    ofstream file("data6.txt");

    file << t << ",";
    file << scalar(V_0[0], V_0[1], V_0[2]) << "\n";
    //file << V_0[0] << ",";
    //file << V_0[1] << ",";
    //file << V_0[2] << "\n";

    //obtain delta values for the 6 variables
    delta_values = h_check(h0, V_0);

    //calculate new h
    h_new = h_optimal(delta_values, h0);

    next_step(h_new, V_0);

    t = t + h_new;

    double next_time = 1.0;

    double delta_time = 1.0;

    while (next_time < 11.0){

        //obtain delta values for the 6 variables
        delta_values = h_check(h_new, V_new);

        //calculate new h
        h_new = h_optimal(delta_values, h_new);

        next_step(h_new, V_new);

        t = t + h_new;
        if ( t > next_time) {

            file << t << ",";
            file << scalar(V_new[0], V_new[1], V_new[2]) << "\n";
	          //file << V_new[0] << "," ;
	          //file << V_new[1] << ",";
	          //file << V_new[2] << "\n";

	          next_time = next_time + delta_time;
        }

     }

}

double scalar(double x, double y, double z){
        s = pow( pow(x, 2.) + pow(y, 2.) + pow(z, 2.), 0.5 );

        return s;
}


int main() {

    a = semimajor(Period_days);
    G_dim = (G* pow(T, 2.0) * Mstar_kg) / pow(a, 3.0); //dimensionless gravitational constant

    h0 = 0.001;

    //Define initial position in dimensionless units

    double x0 = 1.0* pow(1-beta, 1./3.);
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