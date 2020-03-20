#include "functions.cpp"

using namespace std;

//file creators

void file_creator_phi(double* Pa) {//textfile of gaussian vs Rb
  ofstream myfile ("Text_files/phi.txt");
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
  ofstream myfile ("Text_files/theta.txt");
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

// void file_creator_t(int NT, int NP, int NR, double*** t) {//textfile of optical depth vs Rb
//   stringstream title;
//   title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=" << NR;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
//       char string[15];
//       // myfile << Ra[i] << ",";
//       for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
//         for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
//           myfile << t[i][j][k] << endl; //d[theta][phi][R]
//         }
//       }
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

void file_creator_t_ana(int NT, int NP, int NR, double*** t_ana) {//textfile of optical depth vs Rb
  stringstream title;
  title << "Text_files/opticaldepth_ana_nop=" << noparticles ;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << t_ana[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

void file_creator_t_num(int NT, int NP, int NR, double*** t_num) {//textfile of optical depth vs Rb
  stringstream title;
  title << "Text_files/opticaldepth_num_nop=" << noparticles;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << t_num[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file";

}

// void file_creator_total_mass(int NT, int NP, int NR, double*** total_mass) {//textfile of optical depth vs Rb
//   stringstream title;
//   title << "/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/mass=" << NR;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
//       char string[15];
//       // myfile << Ra[i] << ",";
//       for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
//         for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
//           myfile << total_mass[i][j][k] << endl; //d[theta][phi][R]
//         }
//       }
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
//
// }

// void file_creator_d(double*** d, int NT, int NP, int NR) {//textfile of density vs Rb
//   stringstream title;
//   title << "Text_files/density_ana_nop=" << noparticles;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
//       char string[15];
//       // myfile << Ra[i] << ",";
//       for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
//         for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
//           myfile << d[i][j][k] << endl; //d[theta][phi][R]
//         }
//       }
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }
//
// void file_creator_den(double*** den, int NT, int NP, int NR) {//textfile of density vs Rb
//   stringstream title;
//   title << "Text_files/density_num_nop=" << noparticles;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
//       char string[15];
//       // myfile << Ra[i] << ",";
//       for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
//         for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
//           myfile << den[i][j][k] << endl; //d[theta][phi][R]
//         }
//       }
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }
//
// void file_creator_T_vec(double* T_vec, int NT, int NP, int NR) {//textfile of density vs Rb
//   stringstream title;
//   title << "Text_files/T_vec_nop=" << noparticles;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (int x=0; x<noparticles; x++) {
//       char string[15];
//       myfile << T_vec[x] << endl;
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }
//
// void file_creator_P_vec(double* P_vec, int NT, int NP, int NR) {
//   stringstream title;
//   title << "Text_files/P_vec_nop=" << noparticles;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (int x=0; x<noparticles; x++) {
//       char string[15];
//       myfile << P_vec[x] << endl;
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }
//
// void file_creator_R_vec(double* R_vec, int NT, int NP, int NR) {
//   stringstream title;
//   title << "Text_files/R_vec_nop=" << noparticles;
//   ofstream myfile (title.str()+".txt");
//   if (myfile.is_open()) {
//     for (int x=0; x<noparticles; x++) {
//       char string[15];
//       myfile << R_vec[x] << endl;
//     }
//     myfile.close();
//   }
//   else cout << "Unable to open file";
// }
