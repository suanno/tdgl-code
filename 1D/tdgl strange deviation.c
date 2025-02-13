/**
 * @file tdgl_solver.c
 * @brief Solves the time-dependent Ginzburg-Landau (TDGL) equation using the Crank-Nicholson integration scheme.
 *
 * This program numerically integrates the TDGL equation in Fourier space using the Crank-Nicholson scheme,
 * with a time-dependent control parameter C(t). It is possible to define C(t) either 
 * - as a sinusoidal function C(t)=Cave+Ampl*sin(2\pi t/T) where Cave, Ampl, T are specified in the CMD
 *   in this case, the time discretization must be also specified in the CMD.
 * - by loading the values of C(t) from a file "fileCin.dat" in the work directory.
 * 	 In this case, at each step dt is calculated as the distance between the consecutive values of the time specified in "fileCin.dat".
 * The code tracks various physical observables over time. Results are saved for post-processing and analysis
 *
 * ### Usage
 * 	The programs requires to specify a "work directory". This is the directory where it is stored:
 * 	- the state of the system "state.dat". At the end of the run of the present code, it is overwritten.
 *    a further run of the code will CONTINUE the simulation, loading "state.dat" as the initial state.
 * 	- a copy of the initial state (at t=0)
 *  - the files containing the values of the observables which are tracked in time. A file for each observable.
 *  BEFORE the first run, you NEED to
 *  - prepare the initial state, running a code in the /initialization/ directory
 *  - prepare the file "fileCin.dat" defining the values of C(t), if you want to specify C(t) like this and not as an oscillation with the parameter specified in the CMD
 *  The program requires at least two command-line arguments:
 *  - `tspan` (double): Duration of the simulation in time units.
 *  - `simul_path` (string): Path to the folder containing the initial state file (`state.dat`).
 * 	those two parameters are sufficient if the values of C(t) are specified in a "fileCin.dat" file
 *  inside the working directory. Otherwise you should specify also the following parameters:
 *  - `Cave` (double): Average value of the control function `C(t)`.
 *  - `Ampl` (double): Amplitude of `C(t)`.
 *  - `T` (double): Period of `C(t)`. If `T < 0`, `C(t)` remains constant C='Cave'.
 *  - `dt` (double): Time step for the simulation.
 *  - `notsmoothC` (int, optional): If `1`, `C(t)` starts from zero at the beginning of the current evolution.
 *	In the case where C(t) is red from "fileCin.dat", the time discretization dt used at each step is the distance
 *  between the next and current time.
 * 
 * ### Output
 * The simulation saves various files in the work directory, including:
 * - **state.dat**: The final state of the system.
 * - **fileq2Aveout.dat, fileellDW.dat, etc.**: Observables as a function of time.
 * 
 * ### Compilation
 * Compile with:
 * ```bash
 * gcc tdgl.c observables.c read_write.c -fopenmp -lfftw3 -lm -lfftw3_omp -O2 -o ../../1D/.bin/tdgl
 * ```
 *
 * ### Example Run
 * Run a simulation with $C(t) = 1 + 0.5\sin(2\pi t/10)$ for a time 'tspan'=100 and 'dt'=0.1
 * ```bash
 * ../../1D/.bin/tdgl 100 "../../1D/.saves/tempo/" 1 0.5 10 0.1
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

#define MAX_BUFFER_SIZE 256
#define D 1 //Dimension of the problem

#define pi 4*atan(1.0)

int main(int argc, char  *argv [ ]){

/*Simulation parameters: The code reads them from previous simulation .dat file*/
int N, i;
double dx, dt, tspan, tmin, tmax;
double* u;

/* Load the state and prepare the C(t)'s array */
char* simul_path;
char* ptr;
int min_args = 2;                   // Minimum Number of required arguments
if (argc <= min_args){
    printf("Not enought arguments");
    return 0;
}
tspan = (double)strtod(argv[1], &ptr);			// Duration of the simulation
simul_path = argv[2];							// Path containing the simulation's folder
/* Load initial state */
char save_dir[MAX_BUFFER_SIZE] = ""; strcat(save_dir, simul_path);
char state_dir[MAX_BUFFER_SIZE] = ""; strcat(state_dir, save_dir); strcat(state_dir, "/state.dat");
N = loadState(state_dir, &u, &dx, &tmin);
/* Define C(t) */
double* C;   //C[i]=C(t_i+dt)
double* t_C; //Time of the C(t) values
double Cprev;
int nloop = 0;
if (argc == min_args + 1){
	/* Read C(t) from file */
	//char save_dir[MAX_BUFFER_SIZE] = "../../2D/.saves/";
	char fileCinName[MAX_BUFFER_SIZE] = ""; strcat(fileCinName, save_dir); strcat(fileCinName, "/fileCin.dat");
	nloop = readC(fileCinName, &C, &t_C, &Cprev, tmin, tspan);
	if (nloop == 0){ //Error (fileCin.dat too short)
		printf("fileCin.dat is TOO SHORT! or Not found!!");
		return 0;
	}
	//printf("\n%d\n",nloop);
}
else{
	/* C(t)=Cbar+Ampl*sin(2pi t/T) with Cbar, Ampl, T AND dt specified in the CMD */
	double Ampl, T, Cave;						// Amplitude, Half of the period and Average value of C(t)
	int notsmoothC = 0;							// 1: C(t) starts from 0 at the beginning of CURRENT/LAST evolution
												// 0: C(t) starts from zero at the beginning of the FISRT EVER evolution
	Cave = (double)strtod(argv[3], &ptr);
	Ampl = (double)strtod(argv[4], &ptr);
	T = (double)strtod(argv[5], &ptr);
	dt = (double)strtod(argv[6], &ptr);
	if (argc > min_args + 1 + 4)
		notsmoothC = atoi(argv[7]);

	nloop = (int)(tspan/dt);
	C = malloc(nloop*sizeof(double));
	t_C = malloc(nloop*sizeof(double));
	/*Define Cprev (C(t=0))*/
	if(T > 0)
		Cprev = Cave + Ampl*sin(2*pi*(tmin*notsmoothC)/T);
	else
		Cprev = Cave;
	/*Define C(t)*/
	for (int loop=0;loop<nloop;loop++){
		t_C[loop] = tmin + (loop + 1)*dt;		/*C[loop] = C(t+dt), so it's NOT C(t)*/
		if(T > 0)
			C[loop]=Cave + Ampl*sin(2*pi*(t_C[loop]-tmin*notsmoothC)/T);
		else							/*T < 0 means you want to keep C = Cave + Ampl costant*/
			C[loop] = Cave;
	}
}

/*Define observables to track in time*/
int num_saves = 1000; /*Save the observable only at num_saves equispaced instants*/
if (nloop < num_saves)
    num_saves = nloop;
int index_saves = 0;
double* Times = malloc(num_saves*sizeof(double)); /*Times of saves*/
double* q2Ave = malloc(num_saves*sizeof(double));
double* ellDW = malloc(num_saves*sizeof(double));
double* structure_fac = malloc(N*sizeof(double));
double* uave = malloc(num_saves*sizeof(double));
double* kink_dist = malloc(num_saves*sizeof(double));
//double* sigma2ave = malloc(num_saves*sizeof(double));
double weight_sum = 0;
int found_kink = 0;

/*State variables*/
double* x = malloc(N*sizeof(double));
double* ufr = malloc(N*sizeof(double));		//Re FFT[u(x)]
double* ufi = malloc(N*sizeof(double));		//Im FFT[u(x)]
double* ufrdt = malloc(N*sizeof(double));		//Re FFT[u(x)]
double* ufidt = malloc(N*sizeof(double));		//Im FFT[u(x)]
double* NL = malloc(N*sizeof(double));
double* NLfr = malloc(N*sizeof(double));
double* NLfi = malloc(N*sizeof(double));
//double* ux = malloc(N*sizeof(double));		//du/dx
//double* uxfr = malloc(N*sizeof(double));
//double* uxfi = malloc(N*sizeof(double));
/*q-space lattice*/
double* ffr = malloc(N*sizeof(double));
double* qfr = malloc(N*sizeof(double));
double* d2coef = malloc(N*sizeof(double));
double* integ_coef = malloc(N*sizeof(double));

/* FFT */
fftw_complex *in, *out;
fftw_plan pf,pb;
   
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
pf = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
pb = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

/*Define q-space lattice*/
for (int i=1; i<(N/2); i++){
ffr[i]=i;
ffr[N-i]=-i;
}
ffr[0]=0;
ffr[N/2]=N/2;	/*N MUST BE EVEN!!!*/
for (int i=0; i<N; i++){
qfr[i]=ffr[i]*2*pi/(dx*N);
}


// coefficients for spectral derivatives (-q^2) array (with already the minus sign!)
for (i = 0; i<N; i++){
d2coef[i]=-qfr[i]*qfr[i];	//-q^2
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


/* EVOLUTION CODE */
tmax = tmin + tspan;
double ttime = tmin;

/*Time loop*/
for (int loop = 0; loop < nloop; loop++){
    dt = t_C[loop]-ttime;
	ttime = t_C[loop]; //so the first time we calculate (and save in the .dat) u(tmin + dt), NOT u(t)

	/*Denominator of Crank-Nicolson*/
	#pragma omp parallel for //seulement pour les grands systèmes
	for (i = 0; i<N; i++){
	integ_coef[i]=1-dt*C[loop]/2-dt*d2coef[i]/2; 	//Note that d2coeff = -q^2 (already contains a minus sign)
											 		//And that C[loop] is C(t+dt), NOT C(t)
													//And the 1/2 is for Crank-Nicolson
	}

	/*Compute FFT of u(x)->U(q)*/
	for(i=0; i<N; i++) {
		in[i][0]=u[i];
		in[i][1]=0.0;
	}
	fftw_execute(pf); // repeat as needed
	for(i=0; i<N; i++) {
		if (abs(ufi[i]-out[i][1])>1e-9){
			printf("Deviation at t=%lf\n", ttime);
			//printf("%lf, %d, %lf; ", ttime, i, ufr[i]-out[i][0]);
			ufr[i]=out[i][0];
			ufi[i]=out[i][1];
		}
	}

	/*Compute FFT of u(x)^3 (u^3 is called NL: Non Linear term)*/
    #pragma omp parallel for //seulement pour les grands systèmes
	for (i=0; i<N; i++){
	NL[i]=u[i]*u[i]*u[i];
	}
    #pragma omp parallel for //seulement pour les grands systèmes
	for(i=0; i<N; i++) {
	in[i][0]=NL[i];
	in[i][1]=0.0;
	}
	fftw_execute(pf); // repeat as needed
    #pragma omp parallel for //seulement pour les grands systèmes
	for(i=0; i<N; i++) {
	NLfr[i]=out[i][0];
	NLfi[i]=out[i][1];
	}

	/*C[loop] = C(t+dt) but we need even C(t)=Cprev*/
	if (loop > 0)
		Cprev = C[loop-1];
	
	/*Crank-Nicolson*/
    #pragma omp parallel for //seulement pour les grands systèmes
	for (i=0; i<N; i++){
		ufr[i]=(ufr[i]*(1+dt*Cprev/2+dt*d2coef[i]/2)-dt*NLfr[i])/integ_coef[i];
		ufi[i]=(ufi[i]*(1+dt*Cprev/2+dt*d2coef[i]/2)-dt*NLfi[i])/integ_coef[i];	
	}

	/*INVERSE FFT*/
    #pragma omp parallel for //seulement pour les grands systèmes
	for(i=0; i<N; i++) {
	in[i][0]=ufr[i];
	in[i][1]=ufi[i];
	}
	fftw_execute(pb); // repeat as needed
    //#pragma omp parallel for //seulement pour les grands systèmes
	for(i=0; i<N; i++) {
	u[i]=out[i][0]/N;
	}

	/*Compute du/dx
	for(i=0; i<N; i++) {
	uxfr[i]=qfr[i]*ufi[i];
	uxfi[i]=-qfr[i]*ufr[i];
	}
	for(i=0; i<N; i++) {
	in[i][0]=uxfr[i];
	in[i][1]=uxfi[i];
	}
	fftw_execute(pb); // repeat as needed
	for(i=0; i<N; i++) {
	ux[i]=out[i][0]/N;
	}
	*/
	
	/*Measure observables
	if (loop >= ((double)nloop/num_saves)*index_saves){
		Times[index_saves] = ttime;	
		q2Ave[index_saves] = calcq2ave(ufr, ufi, d2coef, N, dx);
		ellDW[index_saves] = calcelllDW(ufr, ufi, d2coef, N, dx);
		kink_dist[index_saves] = calckink_dist(x, u, N, dx);
		//sigma2ave[index_saves] = calcaverage_sigma2(x, u, ux, N, dx);
		uave[index_saves] = calcaverage(ufr, N, dx);

		index_saves = index_saves + 1;
	}*/
}
printf("t = %lf\n", ttime);

/* Save "state.dat" file */
writeState(state_dir, u, N, dx, tmax);

/*Save the structure factor of the final state*/
calcstructure_fact(ufr, ufi, N, structure_fac);
/*Save the values taken by C(t) in time.
  They are appendend, so you save its values from t=0
*/
FILE* varfile;

/*Save the measured observables as a function of time*/
save_observable(varfile, save_dir, "fileCout.dat", t_C, C, nloop, 1);
save_observable(varfile, save_dir, "fileq2Aveout.dat", Times, q2Ave, num_saves, 1);
save_observable(varfile, save_dir, "fileellDW.dat", Times, ellDW, num_saves, 1);
save_observable(varfile, save_dir, "fileSq.dat", qfr, structure_fac, N, 0);
save_observable(varfile, save_dir, "filekinkdist.dat", Times, kink_dist, num_saves, 1);
//save_observable(varfile, "filesigma2ave.dat", Times, sigma2ave, num_saves, 1);
save_observable(varfile, save_dir, "fileuave.dat", Times, uave, num_saves, 1);

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
//free(ux);
//free(uxfr);
//free(uxfi);
free(NL);
free(NLfr);
free(NLfi);
free(C);
free(ffr);
free(qfr);
free(d2coef);
free(integ_coef);

free(Times);
free(q2Ave);
free(ellDW);
free(kink_dist);
//free(sigma2ave);
free(uave);

return 0;
}

