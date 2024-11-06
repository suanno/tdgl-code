#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include </opt/intel/composer_xe_2013_sp1.3.174/mkl/include/fftw/fftw3.h>
#include <fftw3.h>
#include <math.h>
#define D 2 //Dimension of the problem
#define MAX_BUFFER_SIZE 256

#define pi  4*atan(1.0)
int main(int argc, char  *argv [ ]){

double lam_min = 0;     /*Cutoff minimum considered lenghtscale*/
double q2_max = (2*pi/lam_min)*(2*pi/lam_min);

double temp;
/*Space Parameters: From previous simulation*/
int N; double dx;
/*Time Parameters: From params.txt*/
double dt, tmin, tmax, tspan;      /*   tmin: Initial time(of the initial state)
                                        tmax: Final time (of the final state)
                                    */
int read_from_top = 0;              /*Read fileCin.dat from t = 0 instead of t = tmin*/
int loop_read = 0;                  /*If fileCin.dat is finished, it continues reading it from top*/

/* Get inputs from the terminal */
char* fileSimul;
char* ptr;
int min_args = 2;     /* Minimum Number of required arguments*/
                      /*tspan, fileSimul*/
if (argc <= min_args){
    printf("Not enought arguments");
    return 0;
}
tspan = (double)strtod(argv[1], &ptr);
fileSimul = argv[2];
char save_dir[MAX_BUFFER_SIZE] = ".saves/";
strcat(save_dir, fileSimul); /* add the extension */
char fileCinName[MAX_BUFFER_SIZE] = "";
strcat(fileCinName, save_dir);
strcat(fileCinName, "/fileCin.dat");
/* Extra arguments */
int extra_args = 1;
if (argc > min_args + extra_args)
    read_from_top = (int)strtod(argv[min_args + extra_args], &ptr);
extra_args = extra_args + 1;
if (argc > min_args + extra_args)
    loop_read = (int)strtod(argv[min_args + extra_args], &ptr);
extra_args = extra_args + 1;



/*  Read Lattice (space) parameters from last simulation
    [They must be the same for consecutive simulations (for consistency) while dt can change] */
FILE *fileinit;
char params_dir[MAX_BUFFER_SIZE] = "";
//char* params_dir = malloc(strlen(save_dir)+1+strlen("params.txt")+1);
//strcpy(params_dir, save_dir);
strcat(params_dir, save_dir);
strcat(params_dir, "/params.txt");

char state_dir[MAX_BUFFER_SIZE] = "";
//char* state_dir = malloc(strlen(save_dir)+1+strlen("state.dat"));
//strcpy(state_dir, save_dir);
strcat(state_dir, save_dir);
strcat(state_dir, "/state.dat");
int seed;

printf("\n%s\n", state_dir);
fileinit = fopen(state_dir, "r");
/*First line contains parameters*/
fscanf(fileinit, "%d %lf %lf\n", &N, &tmin, &dx);
//fclose(fileinit);                 /*NOT close, we will read the state!*/
/*  Read Time parameters (only dt) from params.txt*/
FILE *fileparams;
fileparams = fopen(params_dir, "r");
fscanf(fileparams, "dx = %lf\ndt = %lf\n", &temp, &dt);   /*Read only dt, use dx of the past simulation for Consistency*/
fclose(fileparams);

/* Time parameters */
tmax = tmin + tspan;
double ttime=0;
int nloop = (tspan)/dt;
int loop = 0;

/* Allocation */
int i;
double** h = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		h[i] = malloc(N * sizeof(double));
double** hfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hfr[i] = malloc(N * sizeof(double));
double** hfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hfi[i] = malloc(N * sizeof(double));
double** hdt = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hdt[i] = malloc(N * sizeof(double));
double** hdtfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hdtfr[i] = malloc(N * sizeof(double));
double** hdtfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hdtfi[i] = malloc(N * sizeof(double));
double** hc = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hc[i] = malloc(N * sizeof(double));
double** hcfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hcfr[i] = malloc(N * sizeof(double));
double** hcfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		hcfi[i] = malloc(N * sizeof(double));

double** q2 = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		q2[i] = malloc(N * sizeof(double));
double** integ_coef = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		integ_coef[i] = malloc(N * sizeof(double));
/*q times F{u(x,y)}*/
double** qxhfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		qxhfr[i] = malloc(N * sizeof(double));
double** qxhfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		qxhfi[i] = malloc(N * sizeof(double));
double** qyhfr = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		qyhfr[i] = malloc(N * sizeof(double));
double** qyhfi = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		qyhfi[i] = malloc(N * sizeof(double));
/*Gradient of u*/
double** ghx = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		ghx[i] = malloc(N * sizeof(double));
double** ghy = malloc(N*sizeof(double*));
for(i = 0; i < N; i++)
		ghy[i] = malloc(N * sizeof(double));

double* ffr = malloc(N*sizeof(double));
double* qfr = malloc(N*sizeof(double));
double L = N*dx; // data for (x,y) is [0, L]x[0, L]
//printf("L = %lf\n", L);



double x, y, z;     /*For reading and printing the values of x,y,h*/
//int i;
int j;
int k;

double deca;

/*
int Nst=8;
double savetimes[8]={1,2,3,4,4000,6000,8000,10000};
*/

double time;


/*Read C(t) from file*/
FILE* fileCin;
double* C = malloc(nloop*sizeof(double));
double Cprev;   /*Last value of C in the previous simulation*/
double decatime = -1;
double decainC;

fileCin = fopen(fileCinName, "r");
/*CHECK IF dt is the SAME in params.txt and in fileCin.dat*/
fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
if (fabs(dt - decatime) > dt/1000){
    printf("The time step dt is different in 'fileCin.dat' (%lf) and 'params.txt' (%lf))", decatime, dt);
    return 0;
}
fclose(fileCin);
decatime = -1;


fileCin = fopen(fileCinName, "r");

/*READ C(t) from file*/
/*CONVENTION: C[i]=C(t_i + dt)*/
/*Bring the pointer in the fileCin file to the right row*/
int ntmin = (int)(tmin/dt); /*How many values of C(t) in the file you have to Ignore*/
for (int i = 0; i <= ntmin*(1-read_from_top); i++){
    /*If you reach the end of fileCin [if t_state > t_finalC]*/
    if(fscanf(fileCin, "%lf %lf \n", &decatime, &decainC) == EOF){
		if (loop_read == 1){
            //fclose(fileCin);
        	//fileCin = fopen(fileCinName, "r");
			rewind(fileCin);
            fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
        }else{
            printf("fileCin.dat is too short! Think to enable the loop reading of fileCin.dat");
            return 0;
        }
	}
    Cprev = decainC;
    loop = loop + 1;
}

printf("Reading C(t) from file. From t = %lf\n", decatime);
/*Now read the values of C[i]=C(t_i+dt)*/
for (int i = ntmin*(1-read_from_top); i < ntmin*(1-read_from_top) + nloop; i++){
    /*If you reach the end of fileCin [if t_state > t_finalC]*/
    if(fscanf(fileCin, "%lf %lf \n", &decatime, &decainC) == EOF){
		if (loop_read == 1){
            //fclose(fileCin);
        	//fileCin = fopen(fileCinName, "r");
			rewind(fileCin);
            fscanf(fileCin, "%lf %lf \n", &decatime, &decainC);
        }else{
            printf("fileCin.dat is too short! Think to enable the loop reading of fileCin.dat");
            return 0;
        }
	}
    C[i-ntmin] = decainC;
    loop = loop + 1;
}
printf("Finish reading C(t) from file. Until t = %lf; tmax = %lf\n", decatime, tmax);
fclose(fileCin);


/*Define observables to track in time*/
int num_saves = 1000; /*Save the observable only at num_saves equispaced instants*/
if (nloop < num_saves)
    num_saves = nloop;
int index_saves = 0;
double* Times = malloc(num_saves*sizeof(double)); /*Times of saves*/
//num_saves = nloop;
/*Space average of U*/
FILE *fileAveout;
double* Ave = malloc(num_saves*sizeof(double)); /*Average magnetization*/
double* q2Ave = malloc(num_saves*sizeof(double)); /*Average magnetization*/
/*Radius (of a circular island)*/
FILE *fileRadiout;
double* integ_grad2 = malloc(num_saves*sizeof(double)); /*Average magnetization*/
double* R2 = malloc(num_saves*sizeof(double)); /*Average of R2 weighted on grad2*/
double weight_sum = 0;  /*Sum of the weights*/


/*-----------------------------------------------------------------------*/

FILE *filenonconserve2dh;
FILE *fileq2mean;

char filename[16];

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

/*Load initial state*/
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        fscanf(fileinit, "%lf %lf %lf\n", &x, &y, &z);
        h[i][j]=z;
        /*printf("u[%d][%d] = %.2lf\n", i, j, h[i][j]);*/
    }
}
fclose(fileinit);   /*Only Now: you can eventually close the state file*/

// coefficients for implicit scheme
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        q2[i][j]=- (qfr[i]*qfr[i]+qfr[j]*qfr[j])/(dx*dx);
    }
}


printf("\n Starting evolution...\n");

/* EVOLUTION CODE */
printf("Number time steps going to simulate = %d\n", nloop);
for(loop=0;loop<nloop;loop++) {
    time = tmin + (loop+1)*dt;
    //printf("T = %lf\n", time);
    /*
    printf("\n\nTime: %.2lf\n", time-dt);
    printf("h[%d][%d]*dt = %.2lf\n", i, j, h[i][j]);
    */
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            integ_coef[i][j]=1-dt*C[loop]/2-dt*q2[i][j]/2;
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

    /*C[loop] = C(t+dt) but we need even C(t)=Cprev for Cranck-Nicolson*/
	//if (tmin == 0)
		//Cprev = C[0];   /*There is no precedent value, so we tool C(t)=C(t+dt) for the FISRT STEP (t=0)*/
    if (loop > 0)
		Cprev = C[loop-1];
    /*If loop = 0 and tmin > 0, the value of Cprev is already set to the last of the last simulation*/

    /* Crank-Nicolson */
    #pragma omp parallel for //seulement pour les grands systèmes
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            hdtfr[i][j]=(hfr[i][j]*(1+dt*Cprev/2+dt*q2[i][j]/2)-dt*hcfr[i][j])/integ_coef[i][j];
            hdtfi[i][j]=(hfi[i][j]*(1+dt*Cprev/2+dt*q2[i][j]/2)-dt*hcfi[i][j])/integ_coef[i][j];
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

    //printf("T = %lf\n", time);

    
    /* Measure Observables (instantaneous value) */
    if (loop >= (nloop/num_saves)*index_saves){
        Times[index_saves] = time;
        /*u space average*/
        Ave[index_saves] = 0;
        for(i=0; i<N; i++){
            for(j=0; j<N; j++){
                Ave[index_saves] = Ave[index_saves] + h[i][j];
            }
        }
        Ave[index_saves] = Ave[index_saves]/(N*N);
        /*1) q2 average*/
        q2Ave[index_saves] = 0;
        weight_sum = 0;
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                if (q2[i][j] < q2_max){
                    /*Minus sign because teh variable q2 is the observable -q2*/
                    q2Ave[index_saves] = q2Ave[index_saves] - q2[i][j]*(hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j]);
                    weight_sum = weight_sum + (hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j]);
                    //printf("u[%d][%d] = %.2lf\n", i, j, h[i][j]);
                }
            }
        }
        //We normalize by the dimension
        q2Ave[index_saves] = q2Ave[index_saves]/(weight_sum*D);
        
        /*2) integral of grad2 and integral of R^2*grad2/integral grad2*/
        for (i=0; i<N; i++){
            for (j=0; j<N; j++){
                qxhfr[i][j]= qfr[i]*hfr[i][j]/dx;     /*hdt and not h because we're considering time t+dt*/
                qxhfi[i][j]= qfr[i]*hfi[i][j]/dx;
                qyhfr[i][j]= qfr[j]*hfr[i][j]/dx;
                qyhfi[i][j]= qfr[j]*hfi[i][j]/dx;
            }
        }
        /*INVERSE FFT for (Grad h)_x*/
        #pragma omp parallel for //seulement pour les grands systèmes
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                in[i*N+j][0] = qxhfi[i][j];
                in[i*N+j][1] = -qxhfr[i][j];
            }
        }
        fftw_execute(pb);
        #pragma omp parallel for //seulement pour les grands systèmes
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                ghx[i][j] = out[i*N+j][0]/(N*N);  /*No need of normalizing as the ratio we calculate has Grad both at numerator and denominator*/
            }
        }        
        /*INVERSE FFT for (Grad h)_y*/
        #pragma omp parallel for //seulement pour les grands systèmes
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                in[i*N+j][0] = qyhfi[i][j];
                in[i*N+j][1] = -qyhfr[i][j];
            }
        }
        fftw_execute(pb);
        #pragma omp parallel for //seulement pour les grands systèmes
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                ghy[i][j] = out[i*N+j][0]/(N*N);  /*No need of normalizing ...*/
            }
        }
        /*Now compute the integral of |Grad_x|^2*/
        integ_grad2[index_saves] = 0;
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++) {
                /*We DO NOT take a sqrt for better resolution of the peak position*/
                integ_grad2[index_saves] = integ_grad2[index_saves] + (ghx[i][j]*ghx[i][j] + ghy[i][j]*ghy[i][j])*(dx*dx);
                R2[index_saves] = R2[index_saves] + ((i*dx-L/2)*(i*dx-L/2)+(j*dx-L/2)*(j*dx-L/2))*(ghx[i][j]*ghx[i][j] + ghy[i][j]*ghy[i][j])*(dx*dx);
            }
        }
        //We do not save the integral of the grad2, but the fraction: TotalArea/integGrad2
        R2[index_saves] = R2[index_saves]/integ_grad2[index_saves];
        integ_grad2[index_saves] = (L*L)/integ_grad2[index_saves];

        index_saves = index_saves + 1;
        //printf("%d: %d\n", loop, (nloop/num_saves)*index_saves);
    }
    
}

printf("\n Saving...\n");

/*Save the final state and a slice of it*/
FILE* filefinalstate;


filefinalstate = fopen(state_dir, "w");
filefinalstate = fopen(state_dir, "w");
/*Save parameters N, tmax, dx, dt, seed, Ampl, Thalf*/
fprintf(filefinalstate, "%d %lf %lf\n", N, tmax, dx);
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        x = i*dx;
        y = j*dx;
        z = h[i][j];
        fprintf(filefinalstate, "%.20f %.20f %.20f\n", x, y, z);
    }
}
fclose(filefinalstate);

/*Save a slice of the final state and of its Gradient squared*/
FILE* fileslice;
char slicestate_dir[MAX_BUFFER_SIZE] = "";
strcat(slicestate_dir, save_dir);
strcat(slicestate_dir, "/sliceState");
int len = snprintf(NULL, 0, "%f", tmax);
char *strtmax = malloc(len + 1);
snprintf(strtmax, len + 1, "%f", tmax);
strcat(slicestate_dir, strtmax);
strcat(slicestate_dir, ".dat");
char slicegrad2_dir[MAX_BUFFER_SIZE] = "";
strcat(slicegrad2_dir, save_dir);
strcat(slicegrad2_dir, "/sliceGrad2");
snprintf(strtmax, len + 1, "%f", tmax);
strcat(slicegrad2_dir, strtmax);
strcat(slicegrad2_dir, ".dat");
//Slice of the state
fileslice = fopen(slicestate_dir, "w");
fprintf(fileslice, "%d %lf %lf\n", N, tmax, dx);
j = (int)(N/2);
for (i=0; i<N; i++){
    x = i*dx;
    z = h[i][j];
    fprintf(fileslice, "%.20f %.20f\n", x, z);
}
fclose(fileslice);
//Slice of Grad2
fileslice = fopen(slicegrad2_dir, "w");
fprintf(fileslice, "%d %lf %lf\n", N, tmax, dx);
j = (int)(N/2);
for (i=0; i<N; i++){
    x = i*dx;
    z = (ghx[i][j]*ghx[i][j] + ghy[i][j]*ghy[i][j]);
    fprintf(fileslice, "%.20f %.20f\n", x, z);
}
fclose(fileslice);

//fclose(filefinalstate);

/*Save the values taken by C(t) in time.
  They are appendend, so you save its values from t=0
*/
FILE* fileCout;
double decaout;
char cout_dir[MAX_BUFFER_SIZE] = "";
//char* obs_dir = malloc(strlen(save_dir)+1+strlen("fileCout.dat"));
//strcpy(obs_dir, save_dir);
strcat(cout_dir, save_dir);
strcat(cout_dir, "/fileCout.dat");
fileCout = fopen(cout_dir, "a");
for (loop=0; loop<nloop; loop++){
    ttime = tmin + (loop+1)*dt;
    decaout = C[loop];
    /*printf("C[%d] = %lf\n", loop, C[loop]);*/
    fprintf(fileCout, "%.5f %.20f\n", ttime, decaout);
}
fclose(fileCout);
//free(obs_dir);
/*Save values of u space average in different times
fileAveout = fopen("fileAveout.dat", "a");
for (i=0; i<num_saves; i++){
    time = Times[i];
    decaout = Ave[i];
    fprintf(fileAveout, "%.5f %.20f\n", time, decaout);
}
fclose(fileAveout);
*/
/*Save values of q2 average in different times*/
//obs_dir = malloc(strlen(save_dir)+1+strlen("fileQ2.dat"));
//strcpy(obs_dir, save_dir);
char q2_dir[MAX_BUFFER_SIZE] = "";
strcat(q2_dir, save_dir);
strcat(q2_dir, "/fileQ2.dat");
fileAveout = fopen(q2_dir, "a");
for (i=0; i<num_saves; i++){
    time = Times[i];
    decaout = q2Ave[i];
    fprintf(fileAveout, "%.5f %.20f\n", time, decaout);
}
//free(obs_dir);
fclose(fileAveout);

/*Save values of integral of Grad2 in different times */
//obs_dir = malloc(strlen(save_dir)+1+strlen("fileGrad2.dat"));
char grad2_dir[MAX_BUFFER_SIZE] = "";
strcat(grad2_dir, save_dir);
//strcpy(obs_dir, save_dir);
strcat(grad2_dir, "/fileGrad2.dat");
fileAveout = fopen(grad2_dir, "a");
for (i=0; i<num_saves; i++){
    time = Times[i];
    decaout = integ_grad2[i];
    fprintf(fileAveout, "%.5f %.20f\n", time, decaout);
}
//free(obs_dir);
fclose(fileAveout);

/*Save values of integral of Grad2 in different times */
//obs_dir = malloc(strlen(save_dir)+1+strlen("fileGrad2.dat"));
char R2_dir[MAX_BUFFER_SIZE] = "";
strcat(R2_dir, save_dir);
//strcpy(obs_dir, save_dir);
strcat(R2_dir, "/fileRadiout.dat");
fileAveout = fopen(R2_dir, "a");
for (i=0; i<num_saves; i++){
    time = Times[i];
    decaout = R2[i];
    fprintf(fileAveout, "%.5f %.20f\n", time, decaout);
}
//free(obs_dir);
fclose(fileAveout);

/*Save FFT of the state
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
filefinalstate = fopen("stateFFT.dat", "w");

fprintf(filefinalstate, "%d %lf %lf\n", N, tmax, dx);
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        x = qfr[i]/dx;
        y = qfr[j]/dx;
        z = sqrt(hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j]);
        fprintf(filefinalstate, "%.20f %.20f %.20f\n", x, y, z);
    }
}
fclose(filefinalstate);
*/


//------------------------------------------

fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_free(in);
fftw_free(out);
// clean up threads
fftw_cleanup_threads();

//free(save_dir);
//free(state_dir);
free(C);
free(h);
free(hc);
free(hfr);
free(hfi);
free(hcfr);
free(hcfi);
free(hdt);
free(hdtfr);
free(hdtfi);
free(q2);
free(qfr);
free(ffr);
free(integ_coef);

return 0;
}

