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
/*Generate N unif r.v. with average hmoy and amplitude eps*/


int i;

int N;
double u0;              /*Average value of u(x)*/
double eps, lam;        /*u(x)=u0+eps\sin(2\pi x/lam)*/


char* ptr;
char* simul_path;
int n_args = 5;         /*Number of required arguments*/
                        /*L, u0, simulation name*/
if (argc <= n_args){
    printf("Not enought input arguments");
    return 0;
}
N = (int)strtod(argv[1], &ptr);
u0 = strtod(argv[2], &ptr);
eps = strtod(argv[3], &ptr);
lam = strtod(argv[4], &ptr);
simul_path = argv[5];

double dx = 0.1;
double dt = 0.01;
/*Read parameters from parameters.txt.
double dx, dt, Ampl, Thalf, Cave;
FILE *fileparams;
fileparams = fopen("parameters.dat", "r");
fscanf(fileparams, "dx %lf\ndt %lf\nA %lf\nT %lf\nCave %lf", &dx, &dt, &Ampl, &Thalf, &Cave);
fclose(fileparams);
*/
printf("dx = %lf\n", dx);


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
double* u = malloc(N*sizeof(double));
double x;
//#pragma omp parallel for
for (i=0; i<N; i++){
    u[i] = u0+eps*sin(2*pi*x/lam);
    x = i*dx;
    fprintf(filestate, "%.5f %.20f\n", x, u[i]);
    fprintf(fileinit, "%.5f %.20f\n", x, u[i]);
}
fclose(fileinit);
fclose(filestate);

/*Recreate fileCout of values of C(t) [Progressive
executions of the dynamics will APPEND info]
FILE *file;
file = fopen("fileCout.dat", "w");
fclose(file);
file = fopen("fileAveout.dat", "w");
fclose(file);
file = fopen("fileq2Aveout.dat", "w");
fclose(file);
file = fopen("filegrad2.dat", "w");
fclose(file);
file = fopen("fileellDW.dat", "w");
fclose(file);
file = fopen("fileumax.dat", "w");
fclose(file);
*/


return 0;

}
