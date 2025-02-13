#ifndef HEADER_H
#define HEADER_H
#include <stdio.h>

#define D 2 //Dimension of the problem
#define pi 4*atan(1.0)

double calcq2ave(double** hfr, double** hfi, double** q2, int N);
double calcCauchyCrofton(double** h, int N, double dx);
double calcRadiusCircularIsland(double** h, int N, double dx);
int measureRadiusCircularIsland(double**h, int N, double dx, double*x0, double*u0);
double calcDW(double** ghx, double** ghy, double dx, int N);
double calcellDW(double** hfr, double** hfi, double** q2, int N, double dx);
int calcstructure_fact(double** ufr, double** ufi, int N, double* structure_fac);
//double calcelllDW(double* ufr, double* ufi, double* d2coef, int N, double dx);
//int calcstructure_fact(double* ufr, double* ufi, int N, double* structure_fac);
//double calcaverage_sigma2(double* x, double* u, double* ux, int N, double dx);
//double calcaverage(double* ufr, int N, double dx);
int save_observable(FILE* dest_file, char* save_dir, char* obs_name, double* obsx, double* obsy, int len, int append);
int save_arraylike_observable(FILE* dest_file, char* save_dir, char* obs_name, double* obsx, double** arraylike_obsy, int lenx, int leny, int append);

#endif