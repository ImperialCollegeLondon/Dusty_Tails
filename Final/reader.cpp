#include<iostream>
#include<fstream>
#pragma pack(push)
#pragma pack(1)
#pragma pack(pop)

using namespace std;

struct DataPoint {
   double time;
   long int id;
   double xPos;
   double yPos;
   double zPos;
   double psize;
 }

void read_binary(string file_name){
  file.open(file_name.c_str(),ios::binary);
  if(!file.is_open())
  {
     cerr << "File noT found \n";
     exit(EXIT_FAILURE);
  }
  // ifstream file_to_open(file_name.c_str(), ios::binary);
  // ofstream file_to_write("output.txt");

  DataPoint data;

  file.read((char*) &data,sizeof(data));

  // for(int i = 0; i < 3; i++)
  //     file_to_write.write((char *) &time[i], sizeof(DataPoint));
  //  file_to_write.close();
  // float value;
  //
  //   if(file_to_open.is_open()){
  //     while(!file_to_open.eof()){
  //       file_to_open.read((char*)&value,sizeof(float));
  //       file_to_write << value << "\n";
  //
  //     }
  //   }
}

 int main() {
   read_binary("/Users/annawilson/Downloads/output_optical_depth.bin");
 }
