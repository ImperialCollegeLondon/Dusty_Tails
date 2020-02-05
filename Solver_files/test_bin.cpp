#include <fstream>

using namespace std;

int main(){
ofstream outfile;
outfile.open("test.bin", ios::out | ios::binary);

double time1  = 1.2;
double time2 = 1.8;


long int id1 = 1357867;
long int id2 = 9876543;

double x1 = 1.4;
double x2 = 2.1;
double y1 = 9.8;
double y2 = 2.3;
double z1 = 8.6;
double z2 = 3.7;

outfile.write((char*) &time1, sizeof(double));
outfile.write((char*) &id1, sizeof(long int));
outfile.write((char*) &x1, sizeof(double));
outfile.write((char*) &y1, sizeof(double));
outfile.write((char*) &z1, sizeof(double));

outfile.write((char*) &time2, sizeof(double));
outfile.write((char*) &id2, sizeof(long int));
outfile.write((char*) &x2, sizeof(double));
outfile.write((char*) &y2, sizeof(double));
outfile.write((char*) &z2, sizeof(double));

outfile.close();
}
