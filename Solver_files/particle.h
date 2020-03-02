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
    double p_temp; //particle temperature

};




vector <Particle> add_particles(vector <Particle> particles, long int current_particles,
                   long int total_particles, double total_t);



void solve_particles(double total_t, double end_t, vector <Particle> particles, \
                     long int total_particles, double t_common, double big_step, \
                     long int current_particles);

vector <Particle> rm_particles(vector <Particle> particles);
vector <Particle> rm_particles2(vector <Particle> particles);
