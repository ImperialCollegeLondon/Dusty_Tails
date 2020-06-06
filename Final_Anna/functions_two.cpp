#include "Const_two.h"

using namespace std;

struct DataPoint {
   double time; //in orbits
   long int id; //integer
   double xPos; //in terms of a (dimensionless)
   double yPos; //in terms of a (dimensionless)
   double zPos; //in terms of a (dimensionless)
   double psize; // in cm - of big particle
   double mass; //in grams - of small particle
   void reset() {
        double time = 0.;
        long int id = 0;
        double xPos = 0.;
        double yPos = 0.;
        double zPos = 0.;
        double psize = 0.;
        double mass = 0.;
    }
};

DataPoint data; //data is a struct of type DataPoints
vector<DataPoint> points;
vector<DataPoint> Particle;
map<double, vector<DataPoint> > dataPointsByTime;
map<double, vector<DataPoint> >::iterator it;
map<double, vector<DataPoint> >::iterator dataPointsByTimeInterator;
vector<double> xPositions, yPositions, zPositions, rPositions, pPositions, tPositions, idPoints, timePoints, Mass;
vector<double> xPositionsAtTime, yPositionsAtTime, zPositionsAtTime, rPositionsAtTime, pPositionsAtTime, tPositionsAtTime, idPointsAtTime, MassAtTime, SizeAtTime;
int NumberOfParticles;

//indices
int is(int NT){
  return 0;
}
int ie(int NT){
  return int(NT);
}

int js(int NP){
  return 0;
}
int je(int NP){
  return int(NP);
}

int ks(int NR){
  return 0;
}
int ke(int NR){
  return int(NR);
}

//cell width
void make_gauss(int NT, int NP, int NR) {
  for(i=is(NT); i<=ie(NT-1); i++){
    double interval = (Tmax-Tmin)/(NT-1) * double(i);
    T_x[i] = interval;
    double gauss_new = exp(-pow(T_x[i]-T_mu,2.)/(2.*pow(T_sd,2.)));
    T_g[i] = gauss_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double interval = (Pmax-Pmin)/(NP-1) * double(j);
    P_x[j] = interval;
    double gauss_new = exp(-pow(P_x[j]-P_mu,2.)/(2.*pow(P_sd,2.)));
    P_g[j] = gauss_new;
  }
  for(k=ks(NR); k<=ke(NR-1); k++){
    double interval = (Rmax-Rmin)/(NR-1) * double(k);
    R_x[k] = interval;
    double gauss_new = exp(-pow((R_x[k]-R_mu),2.)/(2.*pow(R_sd,2.)));
    R_g[k] = gauss_new;
  }
}

void find_inv(int NT, int NP, int NR, double T_B, double P_B, double R_B, double* T_g, double* P_g, double* R_g){
  double T_inv_sum = 0.;
  double P_inv_sum = 0.;
  double R_inv_sum = 0.;
  for(i=is(NT); i<=ie(NT-1); i++){
    double inv_new = 1.-(T_B * T_g[i]);
    if (inv_new == 0.){
      cout << "0 in DT vector" << endl; //so that DR isn't 0
    }
    T_inv[i] = inv_new;
    T_inv_sum += inv_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double inv_new = 1.-(P_B * P_g[j]);
    if (inv_new == 0.){
      cout << "0 in DP vector" << endl; //so that DR isn't 0
    }
    P_inv[j] = inv_new;
    P_inv_sum += inv_new;
  }
  for(k=ks(NR); k<=ke(NR-1); k++){
    double inv_new = 1.-(R_B * R_g[k]);
    if (inv_new == 0.){
      cout << "0 in DR vector" << endl; //so that DR isn't 0
    }
    R_inv[k] = inv_new;
    R_inv_sum += inv_new;
  }
  T_suminv = T_inv_sum;
  P_suminv = P_inv_sum;
  R_suminv = R_inv_sum;
}

void find_DR(int NT, int NP, int NR, double T_A, double P_A, double R_A, double* T_inv, double* P_inv, double* R_inv) {
  double DT_sum = 0.;
  double DP_sum = 0.;
  double DR_sum = 0.;
  for(i=is(NT); i<=ie(NT-1); i++){
    double gauss_new = T_A*T_inv[i];
    DT[i] = gauss_new;
    DT_sum += gauss_new;
  }
  for(j=js(NP); j<=je(NP-1); j++){
    double gauss_new = P_A*P_inv[j];
    DP[j] = gauss_new;
    DP_sum += gauss_new;
  }
  for(k=ks(NR); k<=ke(NR-1); k++){
    double gauss_new = R_A*R_inv[k];
    DR[k] = gauss_new;
    DR_sum += gauss_new;
    //cout << DR[k] << endl;
  }
  sumDT = DT_sum;
  sumDP = DP_sum;
  sumDR = DR_sum;
}

//building the grid
void build_grid(int NT, int NP, int NR, double* DT, double* DP, double* DR){
//R grid
  for (i=is(NT)+1; i<=ie(NT); i++){ // defines the cell edges
    Ta_new = Ta[i-1]+DT[i-1];
    Ta[i] = Ta_new;
  }
  for (i=is(NT); i<=ie(NT)-1; i++){ // defines the cell centers
    Tb_new = Ta[i]+DT[i]/2.;
    Tb[i] = Tb_new;
  }
  for (i=is(NT);i<=ie(NT)-1;i++){ // defines the width of a cell
    dTa_new = Ta[i+1]-Ta[i];
    dTa[i] = dTa_new;
  }
  for (i=is(NT)+1;i<=ie(NT)-1;i++){ // defines the width between two cell centers
    dTb_new = Tb[i]-Tb[i-1];
    dTb[i] = dTb_new;
  }

//PHI grid
  for (j=js(NP)+1; j<=je(NP); j++){ // defines the cell edges
    Pa_new = Pa[j-1]+DP[j-1];
    Pa[j] = Pa_new;
  }
  for (j=js(NP); j<=je(NP)-1; j++){ // defines the cell centers
    Pb_new = Pa[j]+DP[j]/2.;
    Pb[j] = Pb_new;
  }
  for (j=js(NP);j<=je(NP)-1;j++){ // defines the width of a cell
    dPa_new = Pa[j+1]-Pa[j];
    dPa[j] = dPa_new;
  }
  for (j=js(NP)+1;j<=je(NP)-1;j++){ // defines the width between two cell centers
    dPb_new = Pb[j]-Pb[j-1];
    dPb[j] = dPb_new;
  }

//THETA grid
  for (k=ks(NR)+1; k<=ke(NR); k++){ // defines the cell edges
    Ra_new = Ra[k-1]+DR[k-1];
    Ra[k] = Ra_new;
  }
  for (k=ks(NR); k<=ke(NR)-1; k++){ // defines the cell centers
    Rb_new = Ra[k]+DR[k]/2.;
    Rb[k] = Rb_new;
  }
  for (k=ks(NR);k<=ke(NR)-1;k++){ // defines the width of a cell
    dRa_new = Ra[k+1]-Ra[k];
    dRa[k] = dRa_new;
  }
  for (k=ks(NR)+1;k<=ke(NR)-1;k++){ // defines the width between two cell centers
    dRb_new = Rb[k]-Rb[k-1];
    dRb[k] = dRb_new;
  }
}

//reading particle positions
void get_positions(string fileName) {
  ifstream is;
  is.open(fileName, ios::binary);

  while(is.good()){
    is.read((char*)&data.time, sizeof(data.time)); //data.time is the time data point of the data struct
    is.read((char*)&data.id, sizeof(data.id));
    is.read((char*)&data.xPos, sizeof(data.xPos));
    is.read((char*)&data.yPos, sizeof(data.yPos));
    is.read((char*)&data.zPos, sizeof(data.zPos));
    is.read((char*)&data.psize, sizeof(data.psize));
    is.read((char*)&data.mass, sizeof(data.mass));

    it = dataPointsByTime.find(data.time);
    dataPointsByTime[data.time].push_back(data);
  }

  is.close();
  for (dataPointsByTimeInterator = dataPointsByTime.begin(); dataPointsByTimeInterator != dataPointsByTime.end(); dataPointsByTimeInterator++) {

    Time = dataPointsByTimeInterator -> first;
    Particle = dataPointsByTimeInterator -> second;
    timePoints.push_back(Time);

    for (int i=0; i < Particle.size(); i++) {
      // xPositions.push_back(Particle[i].xPos);
      // yPositions.push_back(Particle[i].yPos);
      // zPositions.push_back(Particle[i].zPos);
      // idPoints.push_back(Particle[i].id);
      // Mass.push_back(Particle[i].mass*1e25);
      if (TIME - TOLERANCE <= Time && Time < TIME + TOLERANCE) {
        xPositionsAtTime.push_back(Particle[i].xPos);
        yPositionsAtTime.push_back(Particle[i].yPos);
        zPositionsAtTime.push_back(Particle[i].zPos);
        idPointsAtTime.push_back(Particle[i].id);
        MassAtTime.push_back(Particle[i].mass*1e25);
        SizeAtTime.push_back(Particle[i].psize*1e-2);
      }
    }
  }

  NumberOfParticles = int(idPointsAtTime.size());
  for(int i=0; i < NumberOfParticles; i++){
    // rPositions.push_back(sqrt(pow(xPositions[i],2.)+pow(yPositions[i],2.)+pow(zPositions[i],2.)));
    // tPositions.push_back(acos(zPositions[i]/rPositions[i]));
    // pPositions.push_back(atan(yPositions[i]/xPositions[i]));
    rPositionsAtTime.push_back(sqrt(pow(xPositionsAtTime[i],2.)+pow(yPositionsAtTime[i],2.)+pow(zPositionsAtTime[i],2.)));
    tPositionsAtTime.push_back(acos(zPositionsAtTime[i]/rPositionsAtTime[i]));
    pPositionsAtTime.push_back(atan(yPositionsAtTime[i]/xPositionsAtTime[i]));
  }

  // cout << "r min = " << *min_element(rPositionsAtTime.begin(), rPositionsAtTime.end()) << endl;
  // cout << "r max = " << *max_element(rPositionsAtTime.begin(), rPositionsAtTime.end()) << endl;
  // cout << "p min = " << *min_element(pPositionsAtTime.begin(), pPositionsAtTime.end()) << endl;
  // cout << "p max = " << *max_element(pPositionsAtTime.begin(), pPositionsAtTime.end()) << endl;
  // cout << "t min = " << *min_element(tPositionsAtTime.begin(), tPositionsAtTime.end()) << endl;
  // cout << "t max = " << *max_element(tPositionsAtTime.begin(), tPositionsAtTime.end()) << endl;
}

//filling density pointer
double*** density_fill(int NT, int NP, int NR, vector<double> tPositionsAtTime, vector<double> pPositionsAtTime, vector<double> rPositionsAtTime, vector<double> MassAtTime, vector<double> SizeAtTime, double* Ta, double* Pa, double* Ra)
{
  for(int x=0; x<=NumberOfParticles; x++){//iterating through particles
    for(i=is(NT)+1; i<=ie(NT); i++){
      if((Ta[i]<=tPositionsAtTime[x]) && (tPositionsAtTime[x]<Ta[i+1])){
        for(j=js(NP)+1; j<=je(NP); j++){
          if((Pa[j]<=pPositionsAtTime[x]) && (pPositionsAtTime[x]<Pa[j+1])){
            for(k=ks(NR)+1; k<=ke(NR); k++){
              if((Ra[k]<=rPositionsAtTime[x]) && (rPositionsAtTime[x]<Ra[k+1])){
                double den_new = den[i][j][k]+((MassAtTime[x]/(Mstar_kg*1000))/abs(((pow(Ra[k+1],3.)-pow(Ra[k],3.))/3.)*(-cos(Ta[i+1])+cos(Ta[i]))*(dPa[j])));
                den[i][j][k] = den_new;
                double opacity_new = (3./4.)*(1./density_bulk)*(1./0.5e-6)*(Mstar_kg/pow(a,2.));
                // double opacity_new = (3./4.)*(1./density_bulk)*(1./SizeAtTime[x])*(Mstar_kg/pow(a,2.));
                kappa[i][j][k] = opacity_new;
              }
            }
          }
        }
      }
    }
  }
  return den;
}


//filling opacity pointer
// void opacity_fill(int NT, int NP, int NR, vector<double> SizeAtTime){
//   for(i=is(NT)+1; i<=ie(NT); i++){
//     for(j=js(NP)+1; j<=je(NP); j++){
//       for(k=ks(NR)+1; k<=ke(NR); k++){
//         kappa[i][j][k] = (3./4.)*(1./density_bulk)*(1./SizeAtTime[i])*(Mstar_kg/pow(a,2.));
//       }
//     }
//   }
// }

//calculate optical depth
void calculate_optical_depth(int NT, int NP, int NR, double*** kappa, double*** den){
  t_max = 0.;
  phi_max = 0.;
  theta_max = 0.;
  for(i=is(NT)+1; i<=ie(NT); i++){
    for(j=js(NP)+1; j<=je(NP); j++){
      for(k=ks(NR)+1; k<=ke(NR); k++){
        t[i][j][k]=t[i][j][k-1]+(kappa[i][j][k]*den[i][j][k]*dRa[k]);
        // if(t[i][j][k] > t_max){
        //   t_max = t[i][j][k];
        //   phi_max = Pa[j];
        //   theta_max = Ta[i];
        //   r_max = Ra[k];
        // }
        if((Ta[i]<=((M_PI/2)+0.012)) && ((M_PI/2)-0.012 < Ta[i])){
          if((Pa[j]<=(0.012)) && ((-0.012) < Pa[j])){
            t_planet.push_back(t[i][j][k]);
          }
        }
      }
    }
  }
  t_av_sum = accumulate(t_planet.begin(), t_planet.end(), 0.);
  t_av = (t_av_sum)/(t_planet.size());
  cout << TIME << "," << t_av << endl;
  //cout << t_max << "," << phi_max << "," << theta_max << "," << r_max<< endl;
}
