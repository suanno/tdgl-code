#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include </opt/intel/composer_xe_2013_sp1.3.174/mkl/include/fftw/fftw3.h>
#include <fftw3.h>
#include <math.h>
#include "read_write.h"
#include "observables.h"
#define D 2 //Dimension of the problem
#define MAX_BUFFER_SIZE 256

#define pi  4*atan(1.0)
int main(int argc, char  *argv [ ]){
    double tmin, tspan, tmax;

    int i, j, N;
    double dx;
    double** h;

    /* Read parameters from CMD */
    char* fileSimul;                    //Name of the simulation folder (data is red/wrote ONLY there)
    char* ptr;
    int min_args = 1;                   // Minimum Number of required arguments
    if (argc <= min_args){
        printf("Not enought arguments");
        return 0;
    }
    fileSimul = argv[1];

    char save_dir[MAX_BUFFER_SIZE] = ".saves/"; strcat(save_dir, fileSimul); /* add the extension */
    char state_dir[MAX_BUFFER_SIZE] = ""; strcat(state_dir, save_dir); strcat(state_dir, "/state.dat");

    /*Load initial state*/
    N = loadState(state_dir, &h, &dx, &tmin);
    /*State variables*/
    double** hfr = malloc(N*sizeof(double*));
    for(i = 0; i < N; i++)
            hfr[i] = malloc(N * sizeof(double));
    double** hfi = malloc(N*sizeof(double*));
    for(i = 0; i < N; i++)
            hfi[i] = malloc(N * sizeof(double));
    double* ffr = malloc(N*sizeof(double));
    double* qfr = malloc(N*sizeof(double));
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
    /*Observables*/
    double* structure_fac = malloc(N*sizeof(double));
    
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
    /*Save observables*/
    FILE* varfile;
    calcstructure_fact(hfr, hfi, N, structure_fac);
    save_observable(varfile, save_dir, "fileSq.dat", qfr, structure_fac, N, 0);
    
    free(h);
    free(qfr);
    free(ffr);
    free(hfr);
    free(hfi);
    
    return 1;
}