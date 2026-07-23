/**
 * @file tdgl_solver.c
 * @brief Solves the time-dependent Ginzburg-Landau (TDGL) equation using the Crank-Nicholson integration scheme.
 *
 * This program numerically integrates the TDGL equation in Fourier space using the Crank-Nicholson scheme,
 * with a time-dependent control parameter C(t).
 * The values that C(t) assumes MUST be defined in a file "fileCin.dat" that is inside the "work directory".
 * At each step dt is calculated as the distance between the consecutive values of the time specified in "fileCin.dat".
 * The code tracks various physical observables over time. Results are saved for post-processing and analysis
 *
 * ### Usage
 * 	The programs requires to specify a "work directory". This is the directory where it is stored:
 *  - the file defining C(t) "fileCin.dat".
 * 	- the state of the system "state.dat". At the end of the run of the present code, it is overwritten.
 *    a further run of the code will CONTINUE the simulation, loading "state.dat" as the initial state.
 * 	- a copy of the initial state (at t=0)
 *  - the files containing the values of the observables which are tracked in time. A file for each observable.
 *  BEFORE the first run, you NEED to
 *  - prepare the initial state, running a code in the /initialization/ directory
 *  - prepare the file "fileCin.dat" defining the values of C(t).
 *  The program requires two command-line arguments:
 *  - `tspan` (double): Duration of the simulation in time units.
 *  - `simul_path` (string): Path to the folder containing the initial state file (`state.dat`).
 * 
 * ### Output
 * The simulation saves various files in the work directory, including:
 * - **state.dat**: The final state of the system.
 * - **fileq2Aveout.dat, fileellDW.dat, etc.**: Observables as a function of time.
 * 
 * ### Compilation
 * Compile with:
 * ```bash
 * gcc tdgl.c observables.c read_write.c -fopenmp -lfftw3 -lm -lfftw3_omp -O2 -o ../../2D/.bin/tdgl
 * ```
 *
 * ### Example Run
 * Run a simulation in the working directory "../../2D/.saves/tempo/" for a time 'tspan'=100
 * ```bash
 * ../../2D/.bin/tdgl 100 "../../2D/.saves/tempo/"
 * ```
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include </opt/intel/composer_xe_2013_sp1.3.174/mkl/include/fftw/fftw3.h>
#include <fftw3.h>
#include <math.h>
#include "observables.h"
#include "read_write.h"
#define D 2 //Dimension of the problem
#define MAX_BUFFER_SIZE 256

#define pi  4*atan(1.0)
int main(int argc, char  *argv [ ]){

int i, j, k;
double temp;
int N; double dx;                   //Space
double dt, tmin, tmax, tspan;       //Time
                                    /*  tmin: Initial time(of the initial state)
                                        tmax: Final time (of the final state)
                                    */
/* Allocate memory for the state's variables */
double x, y, z;
// h(x) and its FFT
double** h; //Allocated in "readState" function

/* Read parameters from CMD */
char* simul_path;                    //Name of the simulation folder (data is red/wrote ONLY there)
char* ptr;
int min_args = 2;                   // Minimum Number of required arguments
if (argc <= min_args){
    printf("Not enought arguments");
    return 0;
}
tspan = (double)strtod(argv[1], &ptr);
simul_path = argv[2];

//char save_dir[MAX_BUFFER_SIZE] = "../../2D/.saves/";
char save_dir[MAX_BUFFER_SIZE] = ""; strcat(save_dir, simul_path); /* add the extension */

/*Load initial state*/
char state_dir[MAX_BUFFER_SIZE] = ""; strcat(state_dir, save_dir); strcat(state_dir, "/state.dat");
N = loadState(state_dir, &h, &dx, &tmin);

/*Read C(t)*/
char fileCinName[MAX_BUFFER_SIZE] = ""; strcat(fileCinName, save_dir); strcat(fileCinName, "/fileCin.dat");
double* C;   //C[i]=C(t_i+dt)
double* t_C; //Time of the C(t) values
double Cprev;
int nloop = 0;
nloop = readC(fileCinName, &C, &t_C, &Cprev, tmin, tspan);
if (nloop == 0) //Error (fileCin.dat too short)
    return 0;

/*State variables*/
double** hfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hfr[i] = malloc(N * sizeof(double));
double** hfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hfi[i] = malloc(N * sizeof(double));
double** hdt = malloc(N*sizeof(double*));
// h(x) at time t+dt and its FFT
for(i = 0; i < N; i++)
		hdt[i] = malloc(N * sizeof(double));
double** hdtfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hdtfr[i] = malloc(N * sizeof(double));
double** hdtfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hdtfi[i] = malloc(N * sizeof(double));
// h^3(x) and its FFT
double** hc = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hc[i] = malloc(N * sizeof(double));
double** hcfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hcfr[i] = malloc(N * sizeof(double));
double** hcfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hcfi[i] = malloc(N * sizeof(double));
// q^2 array
double** q2 = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		q2[i] = malloc(N * sizeof(double));
double** integ_coef = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		integ_coef[i] = malloc(N * sizeof(double));

// qfr contains the discrete values available for qx (or qy)
double* ffr = malloc(N*sizeof(double));
double* qfr = malloc(N*sizeof(double));



// Size of the (real-space) simulation lattice
double L = N*dx;
/*Build q-space lattice*/
for (i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;
for (i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/N;
}
// Calculate the values of q2
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        q2[i][j]= (qfr[i]*qfr[i]+qfr[j]*qfr[j])/(dx*dx);
    }
}


/*Tracking the observables in time*/
int num_saves = 1000; /*Save the observable only at num_saves equispaced instants*/
if (nloop < num_saves)
    num_saves = nloop;
int index_saves = 0;
num_saves = nloop;

double* Times = malloc(num_saves*sizeof(double)); /*Times of saves*/
double* Cout = malloc(num_saves*sizeof(double)); 
double* q2Ave = malloc(num_saves*sizeof(double)); 
double* totlenght = malloc(num_saves*sizeof(double));
double* totArea = malloc(num_saves*sizeof(double));
double* m2 = malloc(num_saves*sizeof(double));
double* m4 = malloc(num_saves*sizeof(double));


/* FFTW STUFF */
/* Prepare FFT Routines */
// initialize threads for fft
fftw_init_threads();

fftw_complex *in, *out;
fftw_plan pf, pb;
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N);
// create plans with threads
fftw_plan_with_nthreads(omp_get_max_threads());
// plan for forward transform
pf = fftw_plan_dft_2d(N,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
// pf = fftw_plan_dft_r2c_2d(N,N,in,out,FFTW_ESTIMATE);
// plan for backward transform
pb = fftw_plan_dft_2d(N,N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
// pb = fftw_plan_dft_c2r_2d(N,N,in,out,FFTW_ESTIMATE);
//TUNG: end of adding variables for fft


printf("\n Starting evolution...\n");

FILE* observables_file;
/* EVOLUTION CODE */
double time = tmin;
int loop = 0;
printf("Number time steps going to simulate = %d\n", nloop);
for(loop=0;loop<nloop;loop++) {
    dt = t_C[loop]-time;
    //printf("dt = %lf\n",dt);
    time = time + dt;
    #pragma omp parallel for //seulement pour les grands systèmes
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            integ_coef[i][j]=1-dt*C[loop]/2+dt*q2[i][j]/2;
        }
    }

    /*FFT of h(x,y)->H(qx,qy)*/
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            in[i*N+j][0] = h[i][j];
            in[i*N+j][1] = 0.0;
        }
    }
    fftw_execute(pf);
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            hfr[i][j] = out[i*N+j][0];
            hfi[i][j] = out[i*N+j][1];
        }
    }

    /*FFT of h3(x,y)->H3(qx,qy)*/
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            hc[i][j] = h[i][j]*h[i][j]*h[i][j];
            in[i*N+j][0] = hc[i][j];
            in[i*N+j][1] = 0.0;
        }
    }
    fftw_execute(pf);
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            hcfr[i][j] = out[i*N+j][0];
            hcfi[i][j] = out[i*N+j][1];
        }
    }
    
    if (loop > 0)
		Cprev = C[loop-1];
    
    /* Crank-Nicolson */
    #pragma omp parallel for //seulement pour les grands systèmes
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            hdtfr[i][j]=(hfr[i][j]*(1+dt*Cprev/2-dt*q2[i][j]/2)-dt*hcfr[i][j])/integ_coef[i][j];
            hdtfi[i][j]=(hfi[i][j]*(1+dt*Cprev/2-dt*q2[i][j]/2)-dt*hcfi[i][j])/integ_coef[i][j];
        }
    }

    /*INVERSE FFT*/
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            in[i*N+j][0] = hdtfr[i][j];
            in[i*N+j][1] = hdtfi[i][j];
        }
    }
    fftw_execute(pb);
    #pragma omp parallel for //seulement pour les grands systèmes
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            h[i][j] = out[i*N+j][0]/(N*N);
        }
    }

    //printf("%lf \n", time);
    /* Measure Observables (instantaneous value)*/
    if (loop >= ((double)nloop/num_saves)*index_saves){
        Times[index_saves] = time;
        q2Ave[index_saves] = calcq2ave(hfr, hfi, q2, N);
        totlenght[index_saves] = calcCauchyCrofton(h, N, dx);
        totArea[index_saves] = calcTotalArea(h, N, dx);
        m2[index_saves] = calcm2(h, N);
        m4[index_saves] = calcm4(h, N);
        index_saves = index_saves + 1;
    }
    
}

tmax = time;
printf("\n Saving... tmax = %lf\n", tmax);

writeState(state_dir,h,N,dx,tmax);

save_observable(observables_file, save_dir, "fileQ2.dat", Times, q2Ave, num_saves, 1);
save_observable(observables_file, save_dir, "fileTotlenght.dat", Times, totlenght, num_saves, 1);
save_observable(observables_file, save_dir, "fileTotarea.dat", Times, totArea, num_saves, 1);
save_observable(observables_file, save_dir, "filem2.dat", Times, m2, num_saves, 1);
save_observable(observables_file, save_dir, "filem4.dat", Times, m4, num_saves, 1);
//calcstructure_fact(hfr, hfi, N, structure_fac);

fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);
// clean up threads
fftw_cleanup_threads();

free(C);
free(t_C);
free(h);
free(hfr);
free(hfi);
free(hdt);
free(hdtfr);
free(hdtfi);
free(hc);
free(hcfr);
free(hcfi);
free(ffr);
free(qfr);
free(q2);
free(integ_coef);
/*
free(qxhfr);
free(qxhfi);
free(qyhfr);
free(qyhfi);
free(ghx);
free(ghy);

free(q2Ave);
free(radius_island);
free(totlenght);
free(valPhi);
free(valPhi1);
free(valPhi2);
*/
free(Times);

return 0;
}
