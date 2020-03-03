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
vector <double> k1(3), k2(3), k3(3), k4(3), k5(3), k6(3), k7(3);
vector <double> k1d(3), k2d(3), k3d(3), k4d(3), k5d(3), k6d(3), k7d(3);
double ks1, ks2, ks3, ks4, ks5, ks6, ks7;

void k_values(double h, vector <double> V, bool order5, vector <double> &k1, \
    vector <double> &k2, vector <double> &k3, vector <double> &k4, vector <double> &k5, \
    vector <double> &k6, vector <double> &k7, vector <double> &k1d, \
    vector <double> &k2d, vector <double> &k3d, vector <double> &k4d, vector <double> &k5d, \
    vector <double> &k6d, vector <double> &k7d){

    vector <double> Vtemp(3);
    vector <double> Vtempdot(3);
    double s_temp;

    double s0 = V[6];
    vector <double> V0 = {V[0], V[1], V[2]};
    vector <double> V0dot = {V[3], V[4], V[5]};
    //function to obtain several k values of RK-DP method
    //k1 values
    ks1 = h* sublimation(V[6], V[0], V[1], V[2]);


    for (unsigned int i = 0; i < 3; i++){
        k1d[i] = h* acceleration(i, V0[i]- star_pos[i], \
                   V0[i] - planet_pos[i] + r_planet_dim,
                   V);
        k1[i] = h* V0dot[i];

        Vtemp[i] = V0[i] + a21*k1[i];
        Vtempdot[i] = V0dot[i] + a21*k1d[i];
    }
    s_temp = s0 + a21*ks1;
    //k2 values
    ks2 = h* sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);
    for (unsigned int i = 0; i < 3; i++){
        k2d[i] = h* acceleration(i, Vtemp[i] - star_pos[i], \
                   Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                       Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k2[i] = h* Vtempdot[i];

        Vtemp[i] = V0[i] + a31*k1[i] + a32*k2[i];
        Vtempdot[i] = V0dot[i] + a31*k1d[i] + a32*k2d[i];
    }
    s_temp = s0 + a31*ks1 + a32*ks2;

    ks3 = h* sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);
    //k3 values
    for (unsigned int i = 0; i < 3; i++) {
        k3d[i] = h*acceleration(i, Vtemp[i] - star_pos[i], \
                  Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                      Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k3[i] = h* Vtempdot[i];

        Vtemp[i] = V0[i] + a41*k1[i] + a42*k2[i] + a43*k3[i];
        Vtempdot[i] = V0dot[i] + a41*k1d[i] + a42*k2d[i] + a43*k3d[i];

    }
    s_temp = s0 + a41*ks1 + a42*ks2 + a43*ks3;

    ks4 = h* sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);

    //k4 values
    for (unsigned int i = 0; i < 3; i++) {

        k4d[i] = h* acceleration(i, Vtemp[i] - star_pos[i], \
                   Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                       Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k4[i] = h* Vtempdot[i];

        Vtemp[i] = V0[i] + a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i];
        Vtempdot[i]= V0dot[i] + a51*k1d[i] + a52*k2d[i] + a53*k3d[i] + a54*k4d[i];
    }

    s_temp = s0 + a51*ks1 + a52*ks2 + a53*ks3 + a54*ks4;

    ks5 = h * sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);

    for (unsigned int i = 0; i < 3; i++) {
        k5d[i] = h* acceleration(i, Vtemp[i] - star_pos[i], \
            Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k5[i] = h* Vtempdot[i];

        Vtemp[i] = V0[i] + a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i];
        Vtempdot[i] = V0dot[i] + a61*k1d[i] + a62*k2d[i] + a63*k3d[i] + a64*k4d[i] + a65*k5d[i];
        //cout << k5[i] << endl;
    }

    s_temp = s0 + a61*ks1 + a62*ks2 + a63*ks3 + a64*ks4 + a65*ks5;

    ks6 = h * sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);

    for (unsigned int i = 0; i < 3; i++) {
        k6d[i] = h* acceleration(i, Vtemp[i] - star_pos[i], \
                    Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                        Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k6[i] = h* Vtempdot[i];

        Vtemp[i] = V0[i] + a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i];
        Vtempdot[i] = V0dot[i] + a71*k1d[i] + a73*k3d[i] + a74*k4d[i] + a75*k5d[i] + a76*k6d[i];
    }

    s_temp = s0 + a71*ks1 + a73*ks3 + a74*ks4 + a75*ks5 + a76*ks6;

    if (order5 == true) {
      ks7 = h * sublimation(s_temp, Vtemp[0], Vtemp[1], Vtemp[2]);
      for (unsigned int i = 0; i<3; i++) {

        k7d[i] = h* acceleration(i, Vtemp[i] - star_pos[i], \
                    Vtemp[i] - planet_pos[i] + r_planet_dim, {Vtemp[0], Vtemp[1], Vtemp[2], \
                        Vtempdot[0], Vtempdot[1], Vtempdot[2], s_temp});

        k7[i] = h* Vtempdot[i];
    }

    }
}
