
#ifndef READ_WRITE_H
#define READ_WRITE_H
#include <stdio.h>

#define D 1 //Dimension of the problem
#define pi 4*atan(1.0)

//int readstate(char* fileStatePath, double* h);
int readC(char* fileCinPath, double** C, double** t_C, double* Cprev, double tmin, double tspan);
int loadState(char* fileStatePath, double** h, double* dx, double* tmin);
int writeState(char* fileStatePath, double* h, int N, double dx, double tmax);

#endif