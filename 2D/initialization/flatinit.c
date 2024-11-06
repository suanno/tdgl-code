//#include <iostream.h>
//#include <fstream.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>

#define MAX_BUFFER_SIZE 256
#define pi  4*atan(1.0)

int main(int argc, char  *argv [ ]){

int i,j;

double L, dx;
int N;
double u0;              /*Initial flat value*/
char* simul_path;

/* Read CMD parameters */
char *ptr;
int n_args = 3;         /*Number of required arguments*/
                        /*L, u0, simulation name*/
if (argc <= n_args){
    printf("Not enought input arguments");
    return 0;
}
N = (int)strtod(argv[1], &ptr);
u0 = strtod(argv[2], &ptr);
simul_path = argv[3];

dx = 0.1;
/* Read parameters from params.txt
FILE *fileparams;
fileparams = fopen("params.txt", "r");
fscanf(fileparams, "dx = %lf\ndt = %lf", &dx, &dt);
fclose(fileparams); 
*/
printf("dx = %lf\n", dx);

/* Prepare the save folder */
double x, y, r;

//char save_dir[MAX_BUFFER_SIZE] = "../../2D/.saves/";
char save_dir[MAX_BUFFER_SIZE] = ""; strcat(save_dir, simul_path); /* add the extension */
mkdir(save_dir, 0700);
/*Prepare initial state*/
FILE* filestate;
char state_dir[MAX_BUFFER_SIZE] = ""; strcat(state_dir, save_dir); strcat(state_dir, "/state.dat");
filestate = fopen(state_dir, "w");
/*Backup the initial state (init.dat)*/
char init_dir[MAX_BUFFER_SIZE] = ""; strcat(init_dir, save_dir); strcat(init_dir, "/init.dat");
FILE* fileinit;
fileinit = fopen(init_dir, "w");
fprintf(fileinit, "%d %lf %lf\n", N, 0.0, dx);
fprintf(filestate, "%d %lf %lf\n", N, 0.0, dx);
//#pragma omp parallel for  /*I want the x,y to be SORTED in the state.dat file. So no parallel!*/
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        x = i*dx;
        y = j*dx;
        fprintf(fileinit, "%.5f %.5f %.20f\n", x, y, u0);
        fprintf(filestate, "%.5f %.5f %.20f\n", x, y, u0);
    }
}
printf("State prepared at: %s\n", state_dir);
fclose(fileinit);
fclose(filestate);
/* Copy parameters file 
char params_dir[MAX_BUFFER_SIZE] = ""; strcat(params_dir, save_dir); strcat(params_dir, "/params.txt");
fileinit = fopen(params_dir, "w");
fprintf(fileinit, "dx = %lf\ndt = %lf", dx, dt);
*/

/* Prepare observable folders 
char observables[num_obs][20] = {"/fileQ2.dat", "/fileGrad2.dat", "/fileCout.dat", "/fileAveout.dat"};
for (int i = 0; i < num_obs; i++){
    char obs_dir[MAX_BUFFER_SIZE] = ""; strcat(obs_dir, save_dir); strcat(obs_dir, "/params.txt");
    fileinit = fopen(obs_dir, "w");
    fclose(fileinit);
}
*/


return 0;

}
