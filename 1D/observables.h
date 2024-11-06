#ifndef HEADER_H
#define HEADER_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 1
#define pi 4*atan(1.0)
#define BUF_SIZE 65536
#define MAX_BUFFER_SIZE 256

#define D 1 //Dimension of the problem
#define pi 4*atan(1.0)

double calckink_dist(double* x, double* u, int N, double dx);
double calcq2ave(double* ufr, double* ufi, double* d2coef, int N, double dx);
double calcelllDW(double* ufr, double* ufi, double* d2coef, int N, double dx);
int calcstructure_fact(double* ufr, double* ufi, int N, double* structure_fac);
double calcaverage_sigma2(double* x, double* u, double* ux, int N, double dx);
double calcaverage(double* ufr, int N, double dx);
int save_observable(FILE* dest_file, char* save_dir, char* obs_name, double* obsx, double* obsy, int len, int append);

#endif