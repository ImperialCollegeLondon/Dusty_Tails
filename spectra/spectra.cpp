#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cstring>

#include "opacities_s.h"

using namespace std;

#define amu 1.661e-24
#define kb 1.381e-16

#define Rsun_cgs 6.96e+10 //solar radius
#define G_cgs 6.67259e-8 //gravitational constant
#define c_cgs 2.99792458e+10 //speed of light
#define Msun_cgs 1.9885e+33 // solar mass
#define Mearth_cgs 5.972e+27 // earth mass grams
#define Rearth_cgs 6.378e+8 // earth radius cm

#define PI 3.14159
const double gyr = pow(10.,9) * 365. * 24. * 60. * 60.;
double s_0, rho_d;
int h_cells, v_cells;
double Temp , T, mbig,n_mini,m_star, a_p, inclination, r_star;
 
string opacity_dir;
string opac_data;
Opacities opac;
ofstream output_miri("./output_miri.bin", ios::out | ios::binary);
ofstream output_prism("./output_prism.bin", ios::out | ios::binary);

// struct dust_read{
//     double timestamp;
//     long int id;
//     double x_dust;
//     double y_dust;
//     double z_dust;
//     double vx_dust;
//     double vy_dust;
//     double vz_dust;
//     double nmini;
//     double s_dust;
//     double m_dust;
//     double tau_dust;
//     double temp_dust;
//     double kappa_planck;
//     double kappa_dust_abs;
//     double kappa_dust_scat;
    
// };

struct dust {
   double timestamp;
   double phi;
   double x_p;
   double y_p;
   double z_p;
   double x_dp;
   double y_dp;
   double z_dp;
   double m;
   double nmini;
   double kappa;
   double kappa_scat;
   double size;
   double tau;
};

struct dust_read
{
   double timestamp;
   long int id;
   double x_dust;
   double y_dust;
   double z_dust;
   double vx_dust;
   double vy_dust;
   double vz_dust;
   double s_dust;
   double nmini;
   double h_dust;
   double m_dust;
   double temp_dust;
   double tau_dust;
   double kappa_dust_abs;
   double kappa_dust_scat;
   double kappa_planck;
};

template<typename T>
vector<double> linspace(T start_in, T end_in, int num_in){

  vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

void Opacities::read_data(const char *s_lambda_table, const char *s_size_table,
                          const char *s_stellar_abs_table, const char *s_stellar_scat_table,
                          const char *s_particle_abs_table, const char *s_particle_scat_table,
                          bool log_tables_n)
{
  // open files

  FILE *f_lambda_table = fopen(s_lambda_table, "r");
  FILE *f_size_table = fopen(s_size_table, "r");
  FILE *f_stellar_abs_table = fopen(s_stellar_abs_table, "r");
  FILE *f_stellar_scat_table = fopen(s_stellar_scat_table, "r");
  FILE *f_particle_abs_table = fopen(s_particle_abs_table, "r");
  FILE *f_particle_scat_table = fopen(s_particle_scat_table, "r");

  // declare variables
  double x;

  // read wavelength and size tables from files
  lambda_table = {};
  while (fscanf(f_lambda_table, "%lf ", &x) != EOF)
  {
      lambda_table.push_back(x);
  }

  size_table = {};
  while (fscanf(f_size_table, "%lf ", &x) != EOF)
  {
      size_table.push_back(x);
  }

  // keep table dimensions
  lambda_table_length = lambda_table.size();
  size_table_length = size_table.size();
  cout << "Read lambda and size tables " << endl;
  cout << "lambda length " << lambda_table_length << endl;
  // read absorption and scattering coefficients to stellar radiation
  stellar_abs_table = {};
  for (int i = 0; i < size_table_length; i++)
  {

      fscanf(f_stellar_abs_table, "%lf", &x);
      stellar_abs_table.push_back(x);
  }

  stellar_scat_table = {};
  for (int i = 0; i < size_table_length; i++)
  {
      fscanf(f_stellar_scat_table, "%lf", &x);
      stellar_scat_table.push_back(x);
  }
  cout << "read stellar abs and scat tables " << endl;
  // read absorption and scattering coefficients to particle own radiation
  particle_abs_table = {};
  for (int i = 0; i < size_table_length; i++)
  {
      vector<double> opac_line = {};
      for (int j = 0; j < lambda_table_length; j++)
      {

          fscanf(f_particle_abs_table, "%lf", &x);

          opac_line.push_back(x);
      }
      particle_abs_table.push_back(opac_line);
  }

  particle_scat_table = {};
  for (int i = 0; i < size_table_length; i++)
  {
      vector<double> opac_line = {};
      for (int j = 0; j < lambda_table_length; j++)
      {
          fscanf(f_particle_scat_table, "%lf", &x);
          opac_line.push_back(x);
      }
      particle_scat_table.push_back(opac_line);
  }

  // set whether the temperature and size tables are log-spaced
  log_tables = log_tables_n;

  // close files
  fclose(f_lambda_table);
  fclose(f_size_table);
  fclose(f_stellar_abs_table);
  fclose(f_stellar_scat_table);
  fclose(f_particle_abs_table);
  fclose(f_particle_scat_table);
}

vector <dust_read> read_data(){
  std::fstream input;
  input.open("./input.bin", std::fstream::in | std::fstream::binary);
  input.seekg(0, ios::end);
  int size=input.tellg();
  input.seekg(0, ios::beg);
  long int total = size/sizeof(dust);
  vector <dust_read> dust_grains_out;
  for(int i = 0; i < total; i++){
       dust_grains_out.push_back(dust_read());
       input.read((char *) &dust_grains_out[i], sizeof(dust_read));
       if (dust_grains_out[i].timestamp == 0) {
        dust_grains_out.pop_back();
        break;
       }
    }
   input.close();
   cout << "successfully read the data" << endl;
   return dust_grains_out;
}


double Opacities::particle_abs(double A, double L) const
{
    return interpolate_2d(particle_abs_table, A, L);
}

double Opacities::particle_scat(double A, double L) const
{
    return interpolate_2d(particle_scat_table, A, L);
}


double Opacities::interpolate_2d(const vector<vector<double>> &opac_table, double A, double L) const
{
    // find the size interval in which to interpolate
    //cout << "size   " << A << endl;
    //cout << "temp   " << T << endl;
    int x = size_table_length/2; // approx. index
    if (log_tables)
        x = (int)floor((log10(A)-log10(size_table[0])) * (size_table_length-1) / (log10(size_table[size_table_length-1])-log10(size_table[0])));
    else
        x = (int)floor((A-size_table[0]) * (size_table_length-1) / (size_table[size_table_length-1]-size_table[0]));
    if (x < 0)
        x = 0;
    if (x >= size_table_length)
        x = size_table_length-2;

    if (A <= size_table[0])
    {
        x = 0;
    }
    else if (A >= size_table[size_table_length-1])
    {
        x = size_table_length-1;
    }
    else
    {
        while (A > size_table[x+1])
            x++;
        while (A < size_table[x])
            x--;
    }

    // find the temperature interval in which to interpolate
    int y = lambda_table_length/2; // approx. index
    if (log_tables)
        y = (int)floor((log10(L)-log10(lambda_table[0])) * (lambda_table_length-1) / (log10(lambda_table[lambda_table_length-1])-log10(lambda_table[0])));
    else
        y = (int)floor((L-lambda_table[0]) * (lambda_table_length-1) / (lambda_table[lambda_table_length-1]-lambda_table[0]));
    if (y < 0)
        y = 0;
    if (y >= lambda_table_length)
        y = lambda_table_length-2;

    if (L <= lambda_table[0])
    {
        y = 0;
    }
    else if (L >= lambda_table[lambda_table_length-1])
    {
        y = lambda_table_length-1;
    }
    else
    {
        while (L > lambda_table[y+1])
            y++;
        while (L < lambda_table[y])
            y--;
    }

    // handling cases where input is outside the limits of temperature, size tables
    if ((A <= size_table[0] || A >= size_table[size_table_length-1]) && !(L <= lambda_table[0] || L >= lambda_table[lambda_table_length-1]))
    {
        
        return opac_table[x][y] + (L-lambda_table[y]) * (opac_table[x][y+1]-opac_table[x][y]) / (lambda_table[y+1]-lambda_table[y]);
    }
    else if (!(A <= size_table[0] || A >= size_table[size_table_length-1]) && (L <= lambda_table[0] || L >= lambda_table[lambda_table_length-1]))
    {
        //cout << "at exception because temperature is outside table limit" << endl;
        return opac_table[x][y] + (A-size_table[x]) * (opac_table[x+1][y]-opac_table[x][y]) / (size_table[x+1]-size_table[x]);
    }
    else if ((A <= size_table[0] || A >= size_table[size_table_length-1]) && (L <= lambda_table[0] || L >= lambda_table[lambda_table_length-1]))
    {
        return opac_table[x][y];
    }

    // bilinear interpolation
    double fQ11 = opac_table[x][y];
    double fQ12 = opac_table[x+1][y];
    double fQ21 = opac_table[x][y+1];
    double fQ22 = opac_table[x+1][y+1];

    double fxy1 = (lambda_table[y+1] - L)/(lambda_table[y+1] - lambda_table[y])*fQ11 + (L - lambda_table[y])/(lambda_table[y+1] - lambda_table[y])*fQ21;
    double fxy2 = (lambda_table[y+1] - L)/(lambda_table[y+1] - lambda_table[y])*fQ12 + (L - lambda_table[y])/(lambda_table[y+1] - lambda_table[y])*fQ22;

    double kappa = (size_table[x+1] - A)/(size_table[x+1] - size_table[x])*fxy1 + (A - size_table[x])/(size_table[x+1] - size_table[x])*fxy2;

    // a check kappa is positive, just in case something went wrong
    if (kappa <= 0.)
    {
        printf("Opacity error: opacity has non-positive value \n");
        printf("L = %e, A = %e \n", L, A);
        exit(-1);
    }

    return kappa;
}


int main(){

using namespace std;

   string planet, composition;
   string opac_data, line;
   double no_p;
   double t0;
   double mdot_read; //mass loss rate of planet in Earth masses per Gyr
   
  
   int in_c = 0;
   fstream input_dust;
	input_dust.open("./spectra/spectra.in", ios::in);

    if (input_dust.is_open()){
    while ( getline (input_dust,line) )
    {
      if (in_c == 0) {
        planet = line.substr(0,10);
        cout << planet << endl;
      }
      if (in_c == 1) {
        composition = line.substr(0,15);
        cout << composition << endl;
      }
      if (in_c == 2) {
         s_0 = stod(line.substr(0,5)) * 1.0e-4;
         cout << "s_0" << s_0 << endl;
      }
      if (in_c == 3) {
         no_p =  stod(line.substr(0,5));
      }
      if (in_c ==4) {
         mdot_read = stod(line.substr(0,5));
      }
      if (in_c ==5){
        t0 = stod(line.substr(0,5));
        cout << "t0 " << t0 << endl;
      }
      in_c = in_c + 1;
      
    }
    input_dust.close();
  } else {
     cout << "oops" << endl;
  }

   double mdot = mdot_read*Mearth_cgs/gyr;
   if (composition.substr(0,5) == "Al2O3") {
      cout << "Dust is composed of Corundum." << endl;
      opac_data = "corundum_K95";
      double A = 7.74e+4; //clausius claperyon relation
      double Bp = 39.3; //clausius claperyon relation
      rho_d = 4.0; //dust density
      double mu = 101.961;
      double alpha = 0.1;
   } else if (composition.substr(0,6)=="MgSiO3"){
      cout << "Dust is composed of Enstatite." << endl;
      opac_data = "enstatite_J98_J94_D95";
      double A = 6.89e+4;
      double Bp = 37.8;
      rho_d = 3.20;
      double mu = 100.389;
      double alpha = 0.1;

   } else if (composition.substr(0,12)=="Mg08Fe12SiO4") {
  cout << "Dust is composed of Olivine (Mg08,Fe12)" << endl;
  opac_data = "Mg08Fe12SiO4_J94_D95";
  double A = 6.53e+4; 
  double Bp = 34.3;
  rho_d = 3.80;
  double mu = 178.538;
  double alpha = 0.1;

   } else if (composition.substr(0, 4) == "OlSL") {
   cout << "Dust is composed of Sri Lanka Olivine (Mg1.56 Fe0.4 Si0.91 O4)" << endl;
   opac_data = "OlivineSL_Z11";
   double A = 6.53e+4;
   double Bp = 34.1;
   rho_d = 3.3;
   double mu = 149.81;
   double alpha = 0.1;
} else {
   cout << "Composition unknown, stopping.";
   abort();
}
   if (planet.substr(0,9) == "KIC1255_b") {
      Temp = 4550.0;
      T = 15.68*60.*60.; //planetary period in seconds
      mbig = (mdot * T * 0.01) / no_p ; // 0.01 dependent on when particles are being thrown out of planet
      n_mini = (mbig*3.0) / (rho_d*4.0*PI*pow(s_0, 3));
      cout << "n_mini " << n_mini << endl;
      m_star = 0.67; //stellar mass in solar masses
      double Rstar = 0.66; //stellar radius in solar radii
      double b_p = 0.6; //planetary impact parameter
      a_p = pow((G_cgs*m_star*Msun_cgs* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in cgs
      inclination = acos((Rstar*Rsun_cgs*b_p)/a_p); //in radians
      cout << "inclination " << inclination << endl;
      r_star = (Rstar*Rsun_cgs)/a_p; //stellar radius in terms of the semimajor axis
   } else {
      Temp = 3830.0;
      T = 9.146*60.*60.; //planetary period in seconds
      mbig = (mdot * T * 0.01) / no_p ; // 0.01 dependent on when particles are being thrown out of planet
      m_star = 0.60; //stellar mass in solar masses
      double Rstar = 0.58; //stellar radius in solar radii
      double b_p = 0.65;
      a_p = pow((G_cgs*m_star*Msun_cgs* pow(T, 2.0))/ (4.0*pow(PI, 2.0)), 1.0/3.0); //semi major axis in cgs
      inclination = acos((Rstar*Rsun_cgs*b_p)/a_p); //in radians
      r_star = (Rstar*Rsun_cgs)/a_p; //stellar radius in terms of the semimajor axis
      cout << "rstar" << r_star << endl;
   }
   vector <double> period_steps = {};
   vector <double> miri_mrs_lambda = {};
   vector <double> nirspec_prism_lambda = {};
   vector <dust_read> particles_read;
   vector <dust> particles;
   particles_read = read_data();
   int counter = 0;

   vector <double> timestamps;
   vector <double> transit_depths = {};
   vector <double> grid_cell_size = {};
   int T_int { static_cast<int> (Temp)};
   string T_int_s = to_string(T_int);


   miri_mrs_lambda = linspace(5.0e-4, 28.3e-4, 6000); //using R~3000
   nirspec_prism_lambda = linspace(0.6e-4, 5.3e-4, 156); //R~100
   opacity_dir = "./opacs_jankovic/calc_dust_opac/"+opac_data+"/opac_";
   
   opac.read_data((opacity_dir+"wavelength.dat").c_str(), (opacity_dir+"sizes.dat").c_str(),
        (opacity_dir+"planck_abs_tstar"+T_int_s+".dat").c_str(), 
        (opacity_dir+"planck_sca_tstar"+T_int_s+".dat").c_str(),
        (opacity_dir+"mie_abs_0.02.dat").c_str(), (opacity_dir+"mie_sca_0.02.dat").c_str(),
        true);
   cout << "read opacity data successfully " << endl;
    counter = 0;
    for( dust_read& p : particles_read) {
        particles.push_back(dust());
      
        particles[counter].timestamp = p.timestamp;
        double key = p.timestamp;
        if (find(timestamps.begin(), timestamps.end(), key) == timestamps.end()) {
           timestamps.push_back(key);
           cout << "timestamp is " << key << endl;
       }
        particles[counter].phi = 2.0*PI * (p.timestamp - t0);
        particles[counter].x_p = p.x_dust * cos(particles[counter].phi) - p.y_dust * sin(particles[counter].phi);
        particles[counter].y_p = p.x_dust * sin(particles[counter].phi) + p.y_dust * cos(particles[counter].phi);
        particles[counter].z_p = p.z_dust;
        particles[counter].x_dp = particles[counter].x_p * sin(inclination) + particles[counter].z_p * cos(inclination);
        particles[counter].y_dp = particles[counter].y_p;
        particles[counter].z_dp = -particles[counter].x_p * cos(inclination) + particles[counter].z_p * sin(inclination);
        particles[counter].m = p.m_dust;
        particles[counter].kappa = p.kappa_planck;
        particles[counter].kappa_scat = p.kappa_dust_scat;
        particles[counter].size = p.s_dust;
        particles[counter].nmini = p.nmini;
        particles[counter].tau = p.tau_dust;

        counter = counter + 1;
        }
      for (int m=0; m<timestamps.size(); m++){
      //if (timestamps[m]>1.295 && timestamps[m]<1.305) {

      for (int i=0; i<miri_mrs_lambda.size(); i++){
        double abs = 0.0;
        
      for (dust& p : particles) {
        
        if (p.x_dp > 0.0 && p.timestamp==timestamps[m] && ((pow(p.y_dp,2.)+pow(p.z_dp,2.))<pow(r_star,2.0)))
         {
            abs = abs + (opac.particle_abs(p.size, miri_mrs_lambda[i]) * p.m * p.nmini*exp(-1.0*p.tau) / 
            (4.0*pow(p.x_dp*a_p,2.)+pow(p.y_dp*a_p,2.)+pow(p.z_dp*a_p,2.)));
            
        }
        }
      
      output_miri.write((char*) &miri_mrs_lambda[i], sizeof(double));
      output_miri.write((char*) &abs, sizeof(double));
      //cout << "lambda " << miri_mrs_lambda[i] << endl;
      //cout << "abs " << abs << endl;
      


      }

      for (int i = 0; i < nirspec_prism_lambda.size(); i++)
      {
      double abs = 0.0;
      for (dust &p : particles)
      {

        if (p.x_dp > 0.0 && p.timestamp == timestamps[m] && ((pow(p.y_dp,2.)+pow(p.z_dp,2.))<pow(r_star,2.0))
         ) {
            abs = abs + (opac.particle_abs(p.size, nirspec_prism_lambda[i]) * p.m * p.nmini * exp(-1.0 * p.tau) /
                         (4.0 * pow(p.x_dp * a_p, 2.) + pow(p.y_dp * a_p, 2.) + pow(p.z_dp * a_p, 2.)));
        }
      }

      output_prism.write((char *)&nirspec_prism_lambda[i], sizeof(double));
      output_prism.write((char *)&abs, sizeof(double));
      // cout << "lambda " << miri_mrs_lambda[i] << endl;
      // cout << "abs " << abs << endl;
      }
      //}
      }
      
          
   return 0;
}

