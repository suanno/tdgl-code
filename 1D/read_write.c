#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 1
#define pi 4*atan(1.0)
#define BUF_SIZE 65536
#define MAX_BUFFER_SIZE 256

int readC(char* fileCinPath, double** C, double** t_C, double* Cprev, double tmin, double tspan){
    /*Read C(t) from file*/
    FILE* fileCin;
    double time_fileCin, C_fileCin;
    //Find the size of C(t) array: nloop
    time_fileCin = 0; C_fileCin = 0;
    int nloop = 0; double tfinal = 0;
    fileCin = fopen(fileCinPath, "r");
    while (time_fileCin <= tmin + tspan){
        /*If you reach the end of fileCin [if t_state > t_finalC]*/
        if(fscanf(fileCin, "%lf %lf \n", &time_fileCin, &C_fileCin) == EOF){
            printf("fileCin.dat is too short! Think to enable the loop reading of fileCin.dat");
            fclose(fileCin);
            return 0;
            
        }
        if (tfinal < time_fileCin)
            tfinal = time_fileCin;
        if (time_fileCin > tmin){
            nloop = nloop + 1;
        }
    }
    fclose(fileCin);

    /*CONVENTION: C[i]=C(t_i + dt)*/
    double* C_temp = malloc(nloop*sizeof(double));
    double* t_C_temp = malloc(nloop*sizeof(double)); 
    //printf("%d\n", nloop);
    Cprev = malloc(sizeof(double));
    /*Bring the pointer in the fileCin file to the right row*/
    time_fileCin = 0; C_fileCin = 0;
    fileCin = fopen(fileCinPath, "r");
    while (time_fileCin <= tmin){
        *Cprev = C_fileCin;
        /*If you reach the end of fileCin [if t_state > t_finalC]*/
        if(fscanf(fileCin, "%lf %lf \n", &time_fileCin, &C_fileCin) == EOF){
                printf("fileCin.dat is too short! Think to enable the loop reading of fileCin.dat");
                fclose(fileCin);
                return 0;
        }
    }
    //printf("diff = %f\n", tmin-time_fileCin);
    C_temp[0] = C_fileCin;
    t_C_temp[0] = time_fileCin;
    printf("Reading C(t) from file. From t = %lf\n", time_fileCin);
    for (int i = 1; i < nloop; i++){
        /*If you reach the end of fileCin [if t_state > t_finalC]*/
        if(fscanf(fileCin, "%lf %lf \n", &time_fileCin, &C_fileCin) == EOF){
                printf("fileCin.dat is too short! Think to enable the loop reading of fileCin.dat");
                fclose(fileCin);
                return 0;
        }
        //printf("A");
        C_temp[i] = C_fileCin;
        t_C_temp[i] = time_fileCin;
    }
    printf("Finish reading C(t) from file. Until t = %lf\n", time_fileCin);
    fclose(fileCin);
    *C = C_temp;
    *t_C = t_C_temp;

    return nloop;
}

int loadState(char* fileStatePath, double** hout, double* dx, double* tmin){
    /* Read parameters from the state's file "state.dat" */
    FILE *fileinit;
    double x, z;
    int N,i,j;
    double* h;
    fileinit = fopen(fileStatePath, "r");
    /*First line contains parameters*/
    fscanf(fileinit, "%d %lf %lf\n", &N, tmin, dx);
    h = malloc(N*sizeof(double));
    /*Load initial state*/
    for (i=0; i<N; i++){
        fscanf(fileinit, "%lf %lf\n", &x, &z);
        h[i] = z;
    }
    *hout = h;
    fclose(fileinit);   /*Only Now: you can eventually close the state file*//*Define observables to track over time*/
    
    return N;
}


int writeState(char* fileStatePath, double* h, int N, double dx, double tmax){
    /*Save the final state and a slice of it*/
    FILE* filefinalstate;
    double x, z;
    filefinalstate = fopen(fileStatePath, "w");
    /*Save parameters N, tmax, dx, dt, seed, Ampl, Thalf*/
    fprintf(filefinalstate, "%d %lf %lf\n", N, tmax, dx);
    for (int i=0; i<N; i++){
        x = i*dx;
        z = h[i];
        fprintf(filefinalstate, "%lf %.17g\n", x, z);
        //printf("%lf, ",z);
    }
    fclose(filefinalstate);

    return 1;
}
