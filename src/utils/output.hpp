#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

double output(
int ntotal,
vector<double>& position,
vector<double>& velocity,
vector<double>& density,
vector<double>& energy,
vector<double>& pressure,
int iteration)
{
  int i;
  FILE *file = NULL;
  char buffer [255];
  snprintf(buffer,sizeof(char) * 255, "./data/data_%i.csv",iteration);
  file = fopen(buffer,"w");
  fprintf(file,"position,velocity,density,energy,pressure\n");
  for (i=0;i<ntotal;i++){
     fprintf(file,"%f, %f, %f, %f, %f ",position[i],velocity[i],density[i],energy[i],pressure[i]);
     fprintf(file,"\n");
  }
  fclose(file);
}


#endif // OUTPUT_HPP
