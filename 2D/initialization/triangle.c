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

#include <math.h>

double triangle(double x, double y,
                double x0, double y0, double l)
{
    const double h = sqrt(3.0) * l / 2.0;

    /* Shift coordinates so the center is at the origin */
    x -= x0;
    y -= y0;

    if ( y >= -h/3.0 &&
         y <= 2.0*h/3.0 &&
         y <=  sqrt(3.0)*(x + l/2.0) - h/3.0 &&
         y <= -sqrt(3.0)*(x - l/2.0) - h/3.0 )
        return 1.0;

    return -1.0;
}

void gaussian_convolution(
    double **u,
    double **v,
    int Nx,
    int Ny,
    double dx,
    double sigma)
{
    int R = (int)ceil(3.0*sigma/dx);

    double norm = 1.0/(2.0*M_PI*sigma*sigma);

    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            double sum = 0.0;

            for(int di=-R;di<=R;di++)
            {
                int ii = i+di;

                if(ii<0 || ii>=Nx)
                    continue;

                for(int dj=-R;dj<=R;dj++)
                {
                    int jj = j+dj;

                    if(jj<0 || jj>=Ny)
                        continue;

                    double x = di*dx;
                    double y = dj*dx;

                    double G =
                        norm*exp(-(x*x+y*y)/(2.0*sigma*sigma));

                    sum += G*u[ii][jj];
                }
            }

            v[i][j] = sum*dx*dx;
        }
    }
}

int main(int argc, char  *argv [ ]){

int i,j;

double L, dx;
int N;
double u0, l;              /*Field value is +-u0 and the triangle edge is l*/
char* simul_path;

/* Read CMD parameters */
char *ptr;
int n_args = 4;         /*Number of required arguments*/
                        /*L, u0, simulation name*/
if (argc <= n_args){
    printf("Not enought input arguments");
    return 0;
}
N = (int)strtod(argv[1], &ptr);
u0 = strtod(argv[2], &ptr);
l = strtod(argv[3], &ptr);
simul_path = argv[4];

dx = 0.1;
L=(double)N*dx;
/* Read parameters from params.txt
FILE *fileparams;
fileparams = fopen("params.txt", "r");
fscanf(fileparams, "dx = %lf\ndt = %lf", &dx, &dt);
fclose(fileparams); 
*/
printf("dx = %lf\n", dx);

/* Prepare the save folder */
double x, y, r;
double **u = malloc(N * sizeof(double *));  // Triangle
double **v = malloc(N * sizeof(double *));  // Convoluted triangle

u[0] = malloc(N * N * sizeof(double));
v[0] = malloc(N * N * sizeof(double));

for (i = 1; i < N; i++) {
    u[i] = u[0] + i * N;
    v[i] = v[0] + i * N;
}
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
double x0 = L/2;
double y0 = L/2;
printf("x0: %lf, y0:%lf, l%lf:\n", x0, y0, l);

for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        x = i*dx;
        y = j*dx;
        u[i][j] = u0*triangle(x,y,x0,y0,l);
    }
}
gaussian_convolution(u, v, N, N, dx, 1);
for (i=0; i<N; i++){
    for (j=0; j<N; j++){
        x = i*dx;
        y = j*dx;
        fprintf(fileinit, "%.5f %.5f %.20f\n", x, y, v[i][j]);
        fprintf(filestate, "%.5f %.5f %.20f\n", x, y, v[i][j]);
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
