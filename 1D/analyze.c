// Solve the TDGL equation
// Integration using the Cranck-Nicholson scheme
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "observables.h"
#include "read_write.h"

#define D 1 //Dimension of the problem

#define pi 4*atan(1.0)

int main(int argc, char  *argv [ ]){
/*This code opens the state file "tdgl_results.dat" and calculates some functions of the data*/

int N, i;
double dx;
double time;
double* u;

/* Load the state*/
char* simul_path;
char* ptr;
int min_args = 1;                   // Minimum Number of required arguments
if (argc <= min_args){
    printf("Not enought arguments");
    return 0;
}
simul_path = argv[1];							// Path containing the simulation's folder
/* Load initial state */
char save_dir[MAX_BUFFER_SIZE] = ""; strcat(save_dir, simul_path);
char state_dir[MAX_BUFFER_SIZE] = ""; strcat(state_dir, save_dir); strcat(state_dir, "/state.dat");

N = loadState(state_dir, &u, &dx, &time);
/*State variables*/
double* x = malloc(N*sizeof(double));
double* ufr = malloc(N*sizeof(double));		//Re FFT[u(x)]
double* ufi = malloc(N*sizeof(double));		//Im FFT[u(x)]
/*q-space lattice*/
double* ffr = malloc(N*sizeof(double));
double* qfr = malloc(N*sizeof(double));
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
save_observable(varfile, save_dir, "fileSq.dat", qfr, structure_fac, N, 0);

/*Clear memory*/
fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);

free(x);
free(u);
free(ufr);
free(ufi);
free(structure_fac);
free(ffr);
free(qfr);


return 0;
}

