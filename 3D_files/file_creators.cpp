#include "Non_uniform_grid_3D_pointers.cpp"

using namespace std;

//file creators

void file_creator_phi(double* Pa) {//textfile of gaussian vs Rb
  ofstream myfile ("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/phi.txt");
  if (myfile.is_open()) {
    for (i=0.; i<=NP; i++) {
      char string[15];

      myfile << Pa[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}

void file_creator_theta(double* Ta) {//textfile of gaussian vs Rb
  ofstream myfile ("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/theta.txt");
  if (myfile.is_open()) {
    for (i=0.; i<=NT; i++) {
      char string[15];

      myfile << Ta[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}

// void file_creator_gaussian(vector<double> g, vector <double> x) {//textfile of gaussian vs Rb
//   ofstream myfile ("gaussian_grid.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<g.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << g[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }
//
// void file_creator_inv(vector<double> DR, vector <double> x) {//textfile of non uniform grid DR vs Rb
//   ofstream myfile ("inverse_gaussian_grid.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<DR.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << DR[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }
//
// void file_creator_DR(vector<double> DR, vector <double> Ra) {//textfile of density vs Rb
//   ofstream myfile ("DR.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<DR.size(); i++) {
//       char string[15];
//
//       myfile << Ra[i] << ",";
//       myfile << DR[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

// void file_creator_xRa(vector<double> x, vector <double> Ra) {//textfile of density vs Rb
//   ofstream myfile ("xRa.txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<x.size(); i++) {
//       char string[15];
//
//       myfile << x[i] << ",";
//       myfile << Ra[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

void file_creator_t(int NT, int NP, int NR, double*** t) {//textfile of optical depth vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/opticaldepth_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << t[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}


void file_creator_total_mass(int NT, int NP, int NR, double*** total_mass) {//textfile of optical depth vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/mass=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << total_mass[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_d(double*** d, double* Ta, double* Pa, double* Ra, int NT, int NP, int NR) {//textfile of density vs Rb
  stringstream title;
  title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/density_NR=" << NR;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << d[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}

// void file_creator_den(double*** den, double* Ta, double* Pa, double* Ra, int NT, int NP, int NR) {//textfile of density vs Rb
//   stringstream title;
//   title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/3D_files/density_num_NR=" << NR;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=is(NT); i<ie(NT); i++) { //looping over THETA coordinates
//       char string[15];
//       // myfile << Ra[i] << ",";
//       for (j=js(NP); j<je(NP); j++){ //looping over PHI coordinates
//         for (k=ks(NR); k<ke(NR); k++){ //looping over R coordinates
//           myfile << den[i][j][k] << endl; //d[theta][phi][R]
//         }
//       }
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }

// void file_creator_gauss(vector<double> gauss, vector <double> Ra, int NR) {//textfile of density vs Rb
//   stringstream title;
//   title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NR=" << NR;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=0.; i<gauss.size(); i++) {
//       char string[15];
//
//       myfile << Ra[i] << ",";
//       myfile << gauss[i] << "\n";
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }
