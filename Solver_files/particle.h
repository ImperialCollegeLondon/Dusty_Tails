#include <iostream>
#include <vector>


using namespace std;

class Particle {
   public:
    int id; //unique id of the particle
    vector <double> position; //current position of the particle
    vector <double> velocity; //current velocity of the particle
    double p_size; //size of the particle
    double p_density; //particle density
    double h_updated; //current optimal time step for particle
    //int grid;
    double p_tau; //optical depth

};
