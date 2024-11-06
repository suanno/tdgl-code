// Solve the TDGL equation
// Integration using the Cranck-Nicholson scheme
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "observables.h"
#include "read_write.h"

#define MAX_BUFFER_SIZE 256
#define D 1 //Dimension of the problem

#define pi 4*atan(1.0)

int main(int argc, char  *argv [ ]){

/*Simulation parameters: The code reads them from previous simulation .dat file*/
int N;
double dx, dt, tspan, tmin, tmax;
double x; double* h;

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
N = loadState(state_dir, &h, &dx, &tmin);
/* Define C(t) */
double* C;   //C[i]=C(t_i+dt)
double* t_C; //Time of the C(t) values
double Cprev;
int nloop = 0;
if (argc == min_args){
	/* Read C(t) from file */
	//char save_dir[MAX_BUFFER_SIZE] = "../../2D/.saves/";
	char fileCinName[MAX_BUFFER_SIZE] = ""; strcat(fileCinName, save_dir); strcat(fileCinName, "/fileCin.dat");
	nloop = readC(fileCinName, &C, &t_C, &Cprev, tmin, tspan);
	if (nloop == 0){ //Error (fileCin.dat too short)
		printf("fileCin.dat is TOO SHORT!");
		return 0;
	}
}
else{
	/* C(t)=Cbar+Ampl*sin(2pi t/T) with Cbar, Ampl, T AND dt specified in the CMD */
	double Ampl, T, Cave;						// Amplitude, Half of the period and Average value of C(t)
	Cave = (double)strtod(argv[3], &ptr);
	Ampl = (double)strtod(argv[4], &ptr);
	T = (double)strtod(argv[5], &ptr);
	dt = (double)strtod(argv[6], &ptr);
	nloop = (int)(tspan/dt);
	C = malloc(nloop*sizeof(double));
	

	
}




//C(t) = Cave + Ampl*sin(2pi*t) unless it is required to read C(t) from a file "fileCin.dat"
double Cave, Ampl, Thalf;	// Thalf = T/2; If T < 0 it means C(t) = Cave costant
int notsmoothC = 0;			// 1: C(t) starts from 0 at the beginning of CURRENT/LAST evolution
							// 0: C(t) starts from zero at the beginning of the FISRT EVER evolution




fileinit = fopen("state.dat", "r");
/*First line is for parameters and seed
  tdgl adopts THE SAME parameters (N,dx, dt) that were used in the
  previous evolution with tdglfd or tdgl, unless specified in the command prompt.
*/
fscanf(fileinit, "%d %lf %lf %lf %d %lf %lf %lf\n", &N, &tmin, &dx, &dt, &seed, &Ampl, &Thalf, &Cave);
fclose(fileinit);

/* Get (eventual) new specified parameters from the CMD
 	./tdgl <Deltat> <Ampl> <T> <Cave> <notsmoothC> <dt> 
 	or
 	./tdgl <Deltat> "varfile.dat" 
   Otherwise, simply ./tdgl <Deltat>
*/
char *ptr;
char* fileCname;
int doreadCfromfile = 0;
//printf("argv1 = %lf", atof(argv[1]));
if (argc > 1)
	Deltat = strtod(argv[1], &ptr);
if (argc > 2){
  	Ampl = strtod(argv[2], &ptr);
	if (argc > 3){
		/*Period of the sine*/
		Thalf = strtod(argv[3], &ptr)/2;
		if (argc > 4){
		/*Offset of C(t)*/
		Cave = strtod(argv[4], &ptr);
		}
		if (argc > 5){
			/*Decide wether C(t) must vary SMOOTHLY between
			an evolution and the next one. Or if it shall start from
			sin(0) = 0 in the CURRENT/LAST evolution*/
			if (strtod(argv[5], &ptr) == 1)
				notsmoothC = 1;
		}
		if (argc > 6){
			/*Period of the sine*/
			dt = strtod(argv[6], &ptr);
		}
	}
	else{
		doreadCfromfile = 1;
		fileCname = argv[2];
	}
	//printf("dt = %lf\n", dt);
}

/*The evolutions are consecutive, so the initial time (and state)
 are the finals of the previous evolution
*/
tmax = tmin + Deltat;
double ttime=0;
int nloop=(Deltat)/dt;
int loop;
/*C(t+dt) values*/
double* C = malloc(nloop*sizeof(double));				
double Cprev;

/*Read from file */
int i;
double decainx=0;
double decainu=0;
double decaoutx=0;
double decaoutu=0;
double decaoutC=0;
double decaoutAve=0;
double decainC=0;
double decatime=0;
FILE *varfile;	/*Will use this variable for all files to read/write. So one file at a time is processed.*/

/*Build the function C(t): read from file or C(t)=Cave+Ampl+sin(\pi t/Thalf)*/
/*If the file is shorter than nloop, C(t)
is elonged Periodically*/
if (doreadCfromfile == 1){
	varfile = fopen(fileCname, "r");
	loop = 0;
	while (loop < nloop){
		if(fscanf(varfile, "%lf %lf \n", &decatime, &decainC) == EOF){
			fclose(varfile);
			varfile = fopen(fileCname, "r");
			fscanf(varfile, "%lf %lf \n", &decatime, &decainC);
		}
		C[loop]=decainC;
		//printf("C[%d] = %lf\n", loop, C[loop]);
		loop = loop + 1;
	}
	fclose(varfile);
}else{
	/*Define value of C(t) in time.
		C(t) = Cave + Ampl*sin(pi*t/(T/2))
		NOTE: Time "t" starts at the beginning of the FIRST of
		the serie of consecutive evolutions, and NOT at the beginning
		of the current simulation.
		So the C value at the beginning of the current simulation
		is NOT Ampl*sin(0) = 0 .
		BUT defining C(t) like this, we have the property that C(t)
		is varying SMOOTHLY during the WHOLE (total serie of consecutive evolutions) experiment.
	*/
	/*NOTE: As to calculate u(t+dt) we need C(t+dt), we adopt the following convention
		C[loop] = C(dt*loop + dt)
		So C[0] = C(0 + dt) and so on
	*/
	/*Define Cprev (C(t=0))*/
	if(Thalf > 0)
		Cprev = Cave + Ampl*sin(pi*(tmin*notsmoothC)/Thalf);
	else
		Cprev = Cave;

	for (loop=0;loop<nloop;loop++){
	ttime = tmin + (loop + 1)*dt;		/*C[loop] = C(t+dt), so it's NOT C(t)*/
	if(Thalf > 0){
		C[loop]=Cave + Ampl*sin(pi*(ttime-tmin*notsmoothC)/Thalf);
		}
		else							/*Thalf < 0 means you want to keep C = Cave + Ampl costant*/
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

ttime=0;
loop=0;

char filename[16];


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
fileinit2 = fopen("state.dat", "r");
/*First line is for parameters and seed*/
double Thalf_temp, dt_temp, Cave_temp, Ampl_temp;
fscanf(fileinit2, "%d %lf %lf %lf %d %lf %lf %lf\n", &N, &tmin, &dx, &dt_temp, &seed, &Ampl_temp, &Thalf_temp, &Cave_temp);
/*Now read the initial (smooth) state*/
for (i=0; i<N; i++){
fscanf(fileinit2, "%lf %lf\n", &decainx, &decainu);
x[i]=decainx;
u[i]=decainu;
/*printf("\n\n%lf\n", decainx);*/
}
fclose(fileinit2);

//------------------------------------

// coefficients for spectral derivatives (-q^2) array (with already the minus sign!)
for (i=0; i<N; i++){
d2coef[i]=-qfr[i]*qfr[i];	//-q^2
}

/* EVOLUTION CODE */
for (loop=0; loop < nloop; loop++){
	ttime = tmin + (loop+1)*dt;	/*We calculate u(t+dt) in this loop, 
	so the first time we calculate (and save in the .dat) u(tmin + dt), NOT u(t)*/

	/*Denominator of Crank-Nicolson*/
	for (i=0; i<N; i++){
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
	ufr[i]=out[i][0];
	ufi[i]=out[i][1];
	}

	/*Compute FFT of u(x)^3 (u^3 is called NL: Non Linear term)*/
	for (i=0; i<N; i++){
	NL[i]=u[i]*u[i]*u[i];
	}

	for(i=0; i<N; i++) {
	in[i][0]=NL[i];
	in[i][1]=0.0;
	}
	fftw_execute(pf); // repeat as needed
	for(i=0; i<N; i++) {
	NLfr[i]=out[i][0];
	NLfi[i]=out[i][1];
	}


	/*C[loop] = C(t+dt) but we need even C(t)=Cprev*/
	if (loop > 0)
		Cprev = C[loop-1];
	
	/*Crank-Nicolson*/
	for (i=0; i<N; i++){
		udtfr[i]=(ufr[i]*(1+dt*Cprev/2+dt*d2coef[i]/2)-dt*NLfr[i])/integ_coef[i];
		udtfi[i]=(ufi[i]*(1+dt*Cprev/2+dt*d2coef[i]/2)-dt*NLfi[i])/integ_coef[i];	
	}

	/*INVERSE FFT*/
	for(i=0; i<N; i++) {
	in[i][0]=udtfr[i];
	in[i][1]=udtfi[i];
	}
	fftw_execute(pb); // repeat as needed
	for(i=0; i<N; i++) {
	udt[i]=out[i][0]/N;
	}

	for(i=0; i<N; i++) {
	u[i]=udt[i];
	}

	/*Compute du/dx*/
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
	
	/*Measure observables*/
	if (loop >= ((double)nloop/num_saves)*index_saves){
		Times[index_saves] = ttime;	
		q2Ave[index_saves] = calcq2ave(ufr, ufi, d2coef, N, dx);
		ellDW[index_saves] = calcelllDW(ufr, ufi, d2coef, N, dx);
		kink_dist[index_saves] = calckink_dist(x, u, N, dx);
		//sigma2ave[index_saves] = calcaverage_sigma2(x, u, ux, N, dx);
		uave[index_saves] = calcaverage(ufr, N, dx);

		index_saves = index_saves + 1;
	}
}
printf("t = %lf\n", ttime);

/*Save the final state*/
varfile = fopen("state.dat", "w");
//Save parameters N, tmax, dx, dt, seed, Ampl, Thalf
fprintf(varfile, "%d %.10lf %.10lf %.10lf %d %lf %lf %lf\n", N, tmax, dx, dt, seed, Ampl, Thalf, Cave);
//Save state (x,u) arrays
for (i=0; i<N; i++){
decaoutx=x[i];
decaoutu=u[i];
fprintf(varfile, "%.10f %.20f\n", decaoutx, decaoutu);
}
fclose(varfile);
/*Save the structure factor of the final state*/
calcstructure_fact(ufr, ufi, N, structure_fac);
/*Save the values taken by C(t) in time.
  They are appendend, so you save its values from t=0
*/
varfile = fopen("fileCout.dat", "a");
for (loop=0; loop<nloop; loop++){
ttime = tmin + (loop+1)*dt;
decaoutC = C[loop];
fprintf(varfile, "%.5f %.20f\n", ttime, decaoutC);
}
fclose(varfile);

/*Save the measured observables as a function of time*/
save_observable(varfile, "fileq2Aveout.dat", Times, q2Ave, num_saves, 1);
save_observable(varfile, "fileellDW.dat", Times, ellDW, num_saves, 1);
save_observable(varfile, "fileSq.dat", qfr, structure_fac, N, 0);
save_observable(varfile, "filekinkdist.dat", Times, kink_dist, num_saves, 1);
//save_observable(varfile, "filesigma2ave.dat", Times, sigma2ave, num_saves, 1);
save_observable(varfile, "fileuave.dat", Times, uave, num_saves, 1);

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
free(C);
free(ffr);
free(qfr);
free(d2coef);
free(integ_coef);

free(Times);
free(q2Ave);
free(ellDW);
free(kink_dist);
free(sigma2ave);
free(uave);

return 0;
}

