#include "functions_two.cpp"

using namespace std;

//file creators

void file_creator_phi(double* Pa) {//textfile of gaussian vs Rb
  ofstream myfile ("Text_files_200_035/phi.txt");
  if (myfile.is_open()) {
    for (i=0.; i<=NP; i++) {
      char string[15];

      myfile << Pa[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open phi file";
}

void file_creator_theta(double* Ta) {//textfile of gaussian vs Rb
  ofstream myfile ("Text_files_200_035/theta.txt");
  if (myfile.is_open()) {
    for (i=0.; i<=NT; i++) {
      char string[15];

      myfile << Ta[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open theta file";
}

void file_creator_R(double* Ra) {//textfile of gaussian vs Rb
  ofstream myfile ("Text_files_200_035/R.txt");
  if (myfile.is_open()) {
    for (i=0.; i<=NR; i++) {
      char string[15];

      myfile << Ra[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open R file";
}

void file_creator_rPositions(vector<double> rPositionsAtTime) {
  stringstream title;
  title << "Text_files_200_035/rPositionsAtTime_TIME=" << TIME;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<=rPositionsAtTime.size(); i++) {
      char string[15];

      myfile << rPositionsAtTime[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open rPositionsAtTime file";
}

void file_creator_pPositions(vector<double> pPositionsAtTime) {
  stringstream title;
  title << "Text_files_200_035/pPositionsAtTime_TIME=" << TIME;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<=pPositionsAtTime.size(); i++) {
      char string[15];

      myfile << pPositionsAtTime[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open pPositionsAtTime file";
}

void file_creator_tPositions(vector<double> tPositionsAtTime) {
  stringstream title;
  title << "Text_files_200_035/tPositionsAtTime_TIME=" << TIME;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=0.; i<=tPositionsAtTime.size(); i++) {
      char string[15];

      myfile << tPositionsAtTime[i] << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open tPositionsAtTime file";
}

void file_creator_t(int NT, int NP, int NR, double*** t, double TIME) {//textfile of optical depth vs Rb
  stringstream title;
  title << "Text_files_200_035/optical_depth_TIME=" << TIME;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << t[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open t file";

}

void file_creator_den(double*** den, int NT, int NP, int NR, double TIME) {//textfile of density vs Rb
  stringstream title;
  title << "Text_files_200_035/density_TIME=" << TIME ;
  ofstream myfile (title.str()+".txt");
  if (myfile.is_open()) {
    for (i=is(NT); i<=ie(NT); i++) { //looping over THETA coordinates
      char string[15];
      // myfile << Ra[i] << ",";
      for (j=js(NP); j<=je(NP); j++){ //looping over PHI coordinates
        for (k=ks(NR); k<=ke(NR); k++){ //looping over R coordinates
          myfile << den[i][j][k] << endl; //d[theta][phi][R]
        }
      }
    }
    myfile.close();
  }
  else cout << "Unable to open den file";
}
