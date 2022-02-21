#include <iostream>
#include<iostream>
#include<fstream>
using namespace std;
double size;
struct dust {
   double timestamp;
   long int id;
   double x_dust;
   double y_dust;
   double z_dust;
   double vx_dust;
   double vy_dust;
   double vz_dust;
   double s_dust;
   double h_dust;
   double m_dust;
   double temp_dust;
   double tau_dust;
   double kappa_dust;
};
int main() {
   /*
   ofstream wf("dust.dat", ios::out | ios::binary);
   if(!wf) {
      cout << "Cannot open file!" << endl;
      return 1;
   }
   dust wstu[3];
   wstu[0].id = 1;
   wstu[0].value = 0.1;
   wstu[1].id = 2;
   wstu[1].value = 0.2;
   wstu[2].id = 3;
   wstu[2].value = 0.3;
   for(int i = 0; i < 3; i++)
      wf.write((char *) &wstu[i], sizeof(dust));
   wf.close();
   if(!wf.good()) {
      cout << "Error occurred at writing time!" << endl;
      return 1;
   }
   */
  std::fstream test_file;
  test_file.open("./simulations/KIC1255b_03micro_1mdot_day_025orb_struct_test.bin", std::fstream::in | std::fstream::binary);
   test_file.seekg(0, ios::end);
   size=test_file.tellg();
   test_file.seekg(0, ios::beg);
   cout << "size " << size << endl;
   long int total = size/sizeof(dust);
   dust dust_grains_out[total];
   for(int i = 0; i < total; i++){
      test_file.read((char *) &dust_grains_out[i], sizeof(dust));
      cout << dust_grains_out[i].timestamp << endl;
      cout << dust_grains_out[i].id << endl;
      cout << dust_grains_out[i].x_dust << endl;
      cout << dust_grains_out[i].kappa_dust << endl;
   }
   test_file.close();
   
   /*
   cout<<"Grain details:"<<endl;
   for(int i=0; i < 3; i++) {
      cout << test_s[i].timestamp << endl;
      cout<< test_s[i].id << endl;
      cout << test_s[i].x_dust << endl;
   }
   */
   
   return 0;
}