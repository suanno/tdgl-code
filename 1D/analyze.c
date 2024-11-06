// Solve the TDGL equation
// Integration using the Cranck-Nicholson scheme
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "observables.h"

#define D 1 //Dimension of the problem

#define pi 4*atan(1.0)

int main(int argc, char  *argv [ ]){
/*This code opens the state file "tdgl_results.dat" and calculates some functions of the data*/

int N, i;
double dx;
double tempd; int tempi;

/* Read parameters from "tdgl_results.dat". It contains the final state of the last simulation
											or the initial state, that you just prepared with an
											/initialization/ script*/
FILE *fileinit;
fileinit = fopen("tdgl_result.dat", "r");
/*First line is for parameters and seed
  tdgl adopts THE SAME parameters (N,dx, dt) that were used in the
  previous evolution with tdglfd or tdgl, unless specified in the command prompt.
*/
fscanf(fileinit, "%d %lf %lf %lf %d %lf %lf %lf\n", &N, &tempd, &dx, &tempd, &tempi, &tempd, &tempd, &tempd);
fclose(fileinit);

/*State variables*/
double* x = malloc(N*sizeof(double));
double* u = malloc(N*sizeof(double));
double* ux = malloc(N*sizeof(double));		//du/dx
double* uxfr = malloc(N*sizeof(double));
double* uxfi = malloc(N*sizeof(double));
double* ufr = malloc(N*sizeof(double));		//Re FFT[u(x)]
double* ufi = malloc(N*sizeof(double));		//Im FFT[u(x)]
double* udt = malloc(N*sizeof(double));
double* udtfr = malloc(N*sizeof(double));
double* udtfi = malloc(N*sizeof(double));
double* NL = malloc(N*sizeof(double));
double* NLfr = malloc(N*sizeof(double));
double* NLfi = malloc(N*sizeof(double));
/*q-space lattice*/
double* ffr = malloc(N*sizeof(double));
double* qfr = malloc(N*sizeof(double));
double* d2coef = malloc(N*sizeof(double));
double* integ_coef = malloc(N*sizeof(double));
/*Quantities to calculate*/
double* structure_fac = malloc(N*sizeof(double));

/* FFT */
fftw_complex *in, *out;
fftw_plan pf,pb;
   
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
pf = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
pb = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

/*Define q-space lattice*/
for (i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;	/*N MUST BE EVEN!!!*/
for (i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/(dx*N);
}

/*LOAD INITIAL STATE*/

FILE *fileinit2;
fileinit2 = fopen("tdgl_result.dat", "r");
/*First line is for parameters and seed*/
double Thalf_temp, dt_temp, Cave_temp, Ampl_temp;
double decainx, decainu;
fscanf(fileinit2, "%d %lf %lf %lf %d %lf %lf %lf\n", &N, &tempd, &dx, &tempd, &tempi, &tempd, &tempd, &tempd);
/*Now read the initial (smooth) state*/
for (i=0; i<N; i++){
fscanf(fileinit2, "%lf %lf\n", &decainx, &decainu);
x[i]=decainx;
u[i]=decainu;
/*printf("\n\n%lf\n", decainx);*/
}
fclose(fileinit2);


/*Compute FFT of u(x)->U(q)*/
for(i=0; i<N; i++) {
	in[i][0]=u[i];
	in[i][1]=0.0;
}
fftw_execute(pf); // repeat as needed
for(i=0; i<N; i++) {
	ufr[i]=out[i][0];
	ufi[i]=out[i][1];
}

/*Calculate observables*/
calcstructure_fact(ufr, ufi, N, structure_fac);
/*Save observables*/
FILE* varfile;
save_observable(varfile, "fileSq.dat", qfr, structure_fac, N, 0);

/*Clear memory*/
fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);

free(x);
free(u);
free(ux);
free(ufr);
free(ufi);
free(structure_fac);
free(uxfr);
free(uxfi);
free(udt);
free(udtfr);
free(udtfi);
free(NL);
free(NLfr);
free(NLfi);
free(ffr);
free(qfr);
free(d2coef);
free(integ_coef);


return 0;
}

