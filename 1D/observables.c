#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 1
#define pi 4*atan(1.0)
#define MAX_BUFFER_SIZE 256


int measure_dist(double*u, int N, double dx, double*x0, double*u0){
	/*	Finds 4 around one of the two zeros (the one at x>0)
		it returns them in two arrays of 8 elements x0 and u0.
		In python then you can use this data to do a 3rd degree polynomial interpolation
		
	 */
	double central_plateau,x_next,u_next;
    int i = (int)(N/2);
    central_plateau = u[i];
	while (i < N){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if ((u[i]>0 && central_plateau < 0)||(u[i]<0 && central_plateau > 0)){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
			x0[1] = x_next;
			x0[2] = (i-(int)(N/2))*dx;
			u0[1] = u_next;
			u0[2] = u[i];
			//Take other two points around
			x0[0] = ((i-2)-(int)(N/2))*dx;
			u0[0] = u[i-2];
			x0[3] = ((i+1)-(int)(N/2))*dx;
			u0[3] = u[i+1];
			return 1;
		}
		x_next = (i-(int)(N/2))*dx;
		u_next = u[i];
		i = i + 1;
	}
	return 0;
}

double calckink_dist(double* u, int N, double dx){		
	/* Assuming there are only 2 isolated kinks, it extimates their distance by extimating the position
	   of the zeros of u(x) with a linear fit.
	*/

    double central_plateau,x1,x2,u1,u2,xk,dist;
    int i = (int)(N/2);
	xk = 0; // (right) kink position
    central_plateau = u[i];
	while (i < N){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if ((u[i]>0 && central_plateau < 0)||(u[i]<0 && central_plateau > 0)){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
			x1 = x2;
			x2 = (i-(int)(N/2))*dx;
			u1 = u2;
			u2 = u[i];
			xk = x1 + (u1/(u1-u2))*dx;
			//printf("%lf\n", u1*u2);
			dist = 2*xk;
            return dist;
		}
		x2 = (i-(int)(N/2))*dx;
		u2 = u[i];
		i = i + 1;
	}
    return 0;   //If it does not find any zero, it means that the kinks have a distance comparable with the dx
}

double calcq2ave(double* ufr, double* ufi, double* d2coef, int N, double dx){
	double weight_sum, q2ave;
    weight_sum = 0; q2ave = 0;
	for(int i=0; i<N; i++) {
	    /*Minus sign because the variable q2 is the observable -q2*/
		q2ave = q2ave - d2coef[i]*(ufr[i]*ufr[i] + ufi[i]*ufi[i]);
		weight_sum = weight_sum + (ufr[i]*ufr[i] + ufi[i]*ufi[i]);
	    //printf("u[%d][%d] = %.2lf\n", i, j, h[i][j]);
	}
	q2ave = q2ave/(weight_sum*D);
    return q2ave;
}

double calcInterfaceEnergy(double* ufr, double* ufi, double* d2coef, int N, double dx){
	// Free energy
	double interface_energy, dq;
    interface_energy = 0;
	dq = 2*pi/((double)N*dx);
	for(int i=0; i<N; i++) {
	    /*Minus sign because the variable q2 is the observable -q2*/
		interface_energy = interface_energy - d2coef[i]*(ufr[i]*ufr[i] + ufi[i]*ufi[i])*dq;
	}
	interface_energy = interface_energy/2;
    return interface_energy;
}



double calcIntu2(double* u, int N){
	// Calculate the integral of u(x)^2 normalized by the lenghth of the simulation box
	double sum = 0;
	for(int i=0;i<N;i++) {
		sum = sum + u[i]*u[i];
	}
	return sum/N;
}

double calcIntu(double* u, int N){
	// Calculate the integral of u(x) normalized by the lenghth of the simulation box
	double sum = 0;
	for(int i=0;i<N;i++) {
		sum = sum + u[i];
	}
	return sum/N;
}

// Moments
double calcm2(double* u, int N){
	// Calculate the second moment of phi (not of g!):
	double sum = 0;
	for(int i=0;i<N;i++) {
		sum = sum + u[i]*u[i];
	}
	return 0.5*sum/N;
}

double calcm4(double* u, int N){
	// Calculate the second moment of phi (not of g!):
	double sum = 0;
	for(int i=0;i<N;i++) {
		sum = sum + u[i]*u[i]*u[i]*u[i];
	}
	return sum/N;
}

double calcelllDW(double* ufr, double* ufi, double* d2coef, int N, double dx){
    double grad2 = 0;
    for(int i=0;i<N;i++) {
        grad2 = grad2 - d2coef[i]*(ufr[i]*ufr[i] + ufi[i]*ufi[i])*(2*pi/(N*dx));
	}
	return (N*dx)/grad2;
}

int calcstructure_fact(double* ufr, double* ufi, int N, double* structure_fac){
	for (int i = 0; i < N; i++){
		//NOTICE: You should take an average over many realization to get a smooth curve!!!
		structure_fac[i] = ufr[i]*ufr[i]+ufi[i]*ufi[i];
	}
	return 1;
}

double calcaverage(double* ufr, double* ufi, int N, double dx){
	/*Calculate the average as the q=0 value of the structure factor S(q)*/
	return sqrt(ufr[0]*ufi[0])*sqrt(2*pi)/(N*dx);
}

double calcaverage_sigma2(double* x, double* u, double* ux, int N, double dx){
	/* Extimates the position of each kink by extimating the zeros of u(x) with a linear fit
	   then, assiming the shape of ux^2 is a gaussian, it extimates the sigma2.
	   Returns the average over the kinks.*/

    double xk,x1,x2,uxk,u1,u2,ux1,ux2;
	double sum_sigma2, num_kinks;
    int i = 1;
    u1 = u[0]; x1 = x[0]; ux1 = ux[0];
	u2 = u[i]; x2 = x[i]; ux2 = ux[i];
	num_kinks = 0; sum_sigma2 = 0;
	while (i < N-1){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if (u2*u1 < 0){  //as soon we reach a zero of u(x)
			xk = x1 + (u1/(u1-u2))*dx;
			// Once we have found the kink, let's extimate the sigma of the gaussian approximation
			sum_sigma2 = sum_sigma2 + (x2*x2 - x1*x1 - 2*xk*dx)/(2*log(ux1/ux2));
			num_kinks = num_kinks + 1;
		}
		x1 = x2;
		u1 = u2;
		ux1 = ux2;
		i = i + 1;
		x2 = x[i];
		u2 = u[i];
		ux2 = ux[i];
	}
    return sum_sigma2/num_kinks;
}

double calcnum_kiks(double* u, int N){
	// Measures the number of zeros of u(x)
	int num = 0;
	double u0 = u[0];
	for(int i = 1; i < N; i++){
		if((u[i]<=0 && u0>0)||(u[i]>=0 && u0<0))	//If u(x) changes sign
			num = num + 1;
		u0 = u[i];
	}

	return num;
}

int calcpos_kiks(double* u, int N, double dx, double* x0k){
	// Measure the positions of all kinks and save as time x01 x02 x03 ... x0n
	//where 0,1,2,...,n is the index of the k-th found kink from left to right
	int k = 0;
	double x1,x2,u1,u2;
	u1 = u[0]; x1=0;
	for(int i = 1; i < N; i++){
		u2=u[i]; x2 = i*dx;
		if((u2<=0 && u1>0)||(u2>=0 && u1<0)){	//If u(x) changes sign
			x0k[k] = x1 + (u1/(u1-u2))*dx;
			k = k+1;
		}
		x1=x2; u1=u2;
		i = i + 1;
	}
	// Check if there is a kink between the right and left end of the box (periodic boundary conditions)
	x1 = x2; u1 = u2;
	x2 = N*dx; u2=u[0];
	if((u2<=0 && u1>0)||(u2>=0 && u1<0)){
		x0k[k] = x1 + (u1/(u1-u2))*dx;
		k = k+1;
	}

	return k;	// Number of kinks
}

double calcmin_len(double* u, int N, double dx){
	// Measure the length of the smallest domain
	double len = 0;
	double min_len = (double)N*dx;
	int i_right = N;
	double u0 = u[0];
	// First estimate
	// Border kink
	int i=N-1; u0=u[i];
	while((u[i]<0 && u0<0)||(u[i]>0 && u0>0))	//Right border
		i = i - 1;
	min_len = (N-i)*dx;
	i = 0;
	while((u[i]<0 && u0<0)||(u[i]>0 && u0>0))	//Left border
		i = i + 1;
	min_len = i*dx+min_len;
	//printf("%lf\n",min_len);
	//Other kinks
	len=0; u0=u[i];
	for(int k = i; k < N; k++){
		len = len + dx;
		if((u[k]<=0 && u0>0)||(u[k]>=0 && u0<0)){	//If u(x) changes sign
			if (len < min_len){
				min_len = len;
				i_right = k;
				//printf("a");
			}
			len = dx;
		}
		u0 = u[k];
	}

	// Estimate with better precision, with linear fits
	// Right kink)
	double x1,x2,u1,u2,m,xk_right, xk_left;
	x1 = dx*(i_right-1); x2=dx*(i_right);
	u1 = u[i_right-1]; u2=u[i_right];
	m = (u2-u1)/(x2-x1); xk_right = x1-u1/m;
	// Left kink
	i = i_right-1;
	while (i>=0){
		if ((u[i]>0 && u[i_right-1] < 0)||(u[i]<0 && u[i_right-1] > 0)){	//If u(x) changes sign
			x1=i*dx; u1=u[i];
			x2=(i+1)*dx; u2=u[i+1];
			m = (u2-u1)/(x2-x1); xk_left = x1-u1/m;
			
			//if (min_len<xk_right-xk_left )
			//printf("%lf, %lf\n", min_len, xk_right-xk_left);
			min_len = xk_right-xk_left;
			return min_len;
		}
		i = i - 1;
	}
	//printf("nada\n");
	return min_len;
}

/* For asymmetrical potential */
#define b 2	
double calcPotentialEnergy(double* u, int N, double dx, double a){
	// Potential energy
	double potential_energy;
    potential_energy = 0;
	for(int i=0; i<N; i++) {
		potential_energy = potential_energy + ((1-u[i]*u[i])*(1-u[i]*u[i])*(1-a*tanh(1+b*u[i])))*dx;
	}
    return potential_energy;
}

/* Save */

int save_observable(FILE* dest_file, char* save_dir, char* obs_name, double* obsx, double* obsy, int len, int append){
    double x, y;
    char obs_dir[MAX_BUFFER_SIZE] = "";
    strcat(obs_dir, save_dir);
    strcat(obs_dir, "/");
    strcat(obs_dir, obs_name);
	if (append == 1)
    	dest_file = fopen(obs_dir, "a");
	else
    	dest_file = fopen(obs_dir, "w");
    for (int i=0; i<len; i++){
        x = obsx[i];
        y = obsy[i];
        fprintf(dest_file, "%lf %.17g\n", x, y);
    }
    fclose(dest_file);
    return 1;
}
int save_observable_single_row(FILE* dest_file, char* save_dir, char* obs_name, double obsx, double* obsy, int len, int append){
    double x, y;
    char obs_dir[MAX_BUFFER_SIZE] = "";
    strcat(obs_dir, save_dir);
    strcat(obs_dir, "/");
    strcat(obs_dir, obs_name);
	if (append == 1)
    	dest_file = fopen(obs_dir, "a");
	else
    	dest_file = fopen(obs_dir, "w");
	fprintf(dest_file, "%lf", obsx);
    for (int i=0; i<len; i++){
        fprintf(dest_file, " %.17g", obsy[i]);
    }
	fprintf(dest_file, "\n");
    fclose(dest_file);
    return 1;
}

int save_arraylike_observable(FILE* dest_file, char* save_dir, char* obs_name, double* obsx, double** arraylike_obsy, int lenx, int leny, int append){
    double x;
	double* y;
    char obs_dir[MAX_BUFFER_SIZE] = "";
    strcat(obs_dir, save_dir);
    strcat(obs_dir, "/");
    strcat(obs_dir, obs_name);
	if (append == 1)
    	dest_file = fopen(obs_dir, "a");
	else
    	dest_file = fopen(obs_dir, "w");
    for (int i=0; i<lenx; i++){
        x = obsx[i];
        y = arraylike_obsy[i];
		fprintf(dest_file, "%.5f", x);
		for(int j=0; j < leny; j++)
        	fprintf(dest_file, " %.17g", y[j]);
		fprintf(dest_file, "\n");
    }
    fclose(dest_file);
    return 1;
}