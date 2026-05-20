#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 2
#define pi 4*atan(1.0)
#define BUF_SIZE 65536
#define MAX_BUFFER_SIZE 256

double calcq2ave(double** hfr, double** hfi, double** q2, int N){
    double weight_sum, q2ave;
    int i, j;
    q2ave = 0; weight_sum = 0;
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            /*Minus sign because teh variable q2 is the observable -q2*/
            q2ave = q2ave + q2[i][j]*(hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j]);
            weight_sum = weight_sum + (hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j]);
        }
    }
    q2ave = q2ave/(weight_sum*D);
    return q2ave;
}

int calcTotalArea(double** h, int N, double dx){
    /*Calculate the total area of the u>0 phase by returning the number of pixels with positive u*/
    int num = 0;
    for (int i = 0; i < N; i++){
        for (int j = 1; j < N; j++){
            if (h[j][i] > 0){
                num=num+1;
            }
        }
    }
    return num;
}


double calcPhi(double** h, int N, double dx){
    /*Calculate Phi=sqrt of integral phi^2/L^D*/
    double Phi = 0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Phi = Phi + h[i][j]*h[i][j];
        }
    }
    return sqrt(Phi/(N*N));
}
double calcPhi1(double** h, int N, double dx){
    /*Calculate Phi1=integral phi/L^D*/
    double Phi1 = 0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Phi1 = Phi1 + h[i][j];
        }
    }
    return Phi1/(N*N);
}
double calcPhi2(double** h, int N, double dx){
    /*Calculate Phi2=sqrt of integral (phi-phi_-)/L^D*/
    double Phi2 = 0;
    double phi_ = h[N-1][N-1];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Phi2 = Phi2 + h[i][j];
        }
    }
    return Phi2/(N*N)-phi_;
}

double calcCauchyCrofton(double** u, int N, double dx){
    /*Extimate the (total) interface(s) lenght by using the Cauchy-Croft formula */
    /*N.B.: In the book of "do Carmo, Manfredo" (P.48) there is a pi/4 instead of pi
            the reason is that it takes 4 sets of parallel lines and considers n to be
            the SUM of the number of crossings with ALL the sets of lines.
            Instead, here, we write a formula for the estimate associated to EACH set of lines,
            we sum the estimates and, at the end, we divide by the number of sets (so we get an average)
    */

    int i, j, num_crossings;
    double average_lenght, uprev;
    average_lenght = 0;
    //Horizontal lines, spaced dx
    num_crossings = 0;
    for (i = 0; i < N; i++){
        uprev = u[i][0];
        for (j = 1; j < N; j++){
            if (u[i][j]*uprev < 0){
                num_crossings = num_crossings + 1;
            }
            uprev = u[i][j];
        }
    }
    average_lenght = average_lenght + 0.5*num_crossings*dx*(pi);
    //Vertical lines, spaced dx
    num_crossings = 0;
    for (i = 0; i < N; i++){
        uprev = u[0][i];
        for (j = 1; j < N; j++){
            if (u[j][i]*uprev < 0){
                num_crossings = num_crossings + 1;
            }
            uprev = u[j][i];
        }
    }
    average_lenght = average_lenght + 0.5*num_crossings*dx*(pi);
    //pi/4 lines, spaced sqrt(2)*dx
    num_crossings = 0;
    for (i = 1; i < N; i++){
        uprev = u[0][i];
        j = 1;
        while (j <= i){
            if (u[j][i-j]*uprev < 0){
                num_crossings = num_crossings + 1;
            }
            uprev = u[j][i-j];
            j = j + 1;
        }
    }
    average_lenght = average_lenght + 0.5*num_crossings*(sqrt(2)*dx)*(pi);
    //pi*3/4 lines, spaced sqrt(2)*dx
    num_crossings = 0;
    for (i = 1; i < N; i++){
        uprev = u[N-1][i];
        j = 1;
        while (j <= i){
            if (u[N-1-j][i-j]*uprev < 0){
                num_crossings = num_crossings + 1;
            }
            uprev = u[N-1-j][i-j];
            j = j + 1;
        }
    }
    average_lenght = average_lenght + 0.5*num_crossings*(sqrt(2)*dx)*(pi);

    return average_lenght/4;
}

double calcRadiusCircularIsland(double** h, int N, double dx){
    // Estimate the radius of a circular island (centered at the origin) by estimating the distance between the zeros
    // of u(x.y) along an horizontal line passing through the origin
    double central_plateau,x1,x2,u1,u2,xk;
    double radius = 0;
    int middle_index = (int)(N/2);
    int i = middle_index;
    int found_zero = 0;
	xk = 0; // (right) kink position
    central_plateau = h[i][i];

    //Horizontal line
	while (i < N && found_zero == 0){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if (h[i][middle_index]*central_plateau < 0){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
			x1 = x2;
			x2 = (i-(int)(N/2))*dx;
			u1 = u2;
			u2 = h[i][middle_index];
			xk = x1 + (u1/(u1-u2))*dx;
			//printf("%lf\n", u1*u2);
			radius = radius + xk;
            found_zero = 1;
		}
		x2 = (i-(int)(N/2))*dx;
		u2 = h[i][middle_index];
		i = i + 1;
	}
    //Vertical line
    found_zero = 0;
    i = middle_index;
	while (i < N && found_zero == 0){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if (h[middle_index][i]*central_plateau < 0){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
			x1 = x2;
			x2 = (i-(int)(N/2))*dx;
			u1 = u2;
			u2 = h[middle_index][i];
			xk = x1 + (u1/(u1-u2))*dx;
			//printf("%lf\n", u1*u2);
			radius = radius + xk;
            found_zero = 1;
		}
		x2 = (i-(int)(N/2))*dx;
		u2 = h[middle_index][i];
		i = i + 1;
	}
    return radius/2;   
}

typedef struct {
    int x;
    int y;
} Direction2D;

int calcRadiusCircularIsland_pi4(double** h, int N, double dx, double*list_R, double*list_tantheta){
    // Estimate the radius of a circular island (centered at the origin) by estimating the distance between the zeros
    // of u(x.y) along some line passing through the origin with an angle between zero and pi/4
    double central_plateau,x1,x2,y1,y2,u1,u2,radius,dr,dl;
    int middle_index = (int)(N/2);
    int j,k;
    int found_zero = 0;
    central_plateau = h[middle_index][middle_index];

    //Cycle over the directions
    Direction2D dirs[] = {
        {1,0},
        {6,1},
        {5, 1},
        {3, 1},
        {2, 1},
        {3, 2},
        {5,4},
        {7,4},
        {1,1}
        // Add as many as you like
    };
    int num_pairs = sizeof(dirs) / sizeof(dirs[0]);
    for (int i = 0; i < num_pairs; i++) {
        // Find the radius with a linear interpolation
        found_zero = 0;
        j = middle_index;
        k = middle_index;
        x1=0;
        x2=0;
        while (j < N && k < N && found_zero == 0){
            // Estimate the position of the zero x=xk (y0=0) with a linear fit
            // We store the value of y(=u) until it changes sign, so we can use the previous value to do the Interpolation (linear passing through two points)
            if (h[j][k]*central_plateau < 0){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
                x1 = x2;
                y1 = y2;
                u1 = u2;
                x2 = (j-(int)(N/2))*dx;
                y2 = (k-(int)(N/2))*dx;
                u2 = h[j][k];
                dl = sqrt(dirs[i].x*dirs[i].x+dirs[i].y*dirs[i].y)*dx;
                dr = u1/(u1-u2)*dl;
                radius = sqrt(x1*x1+y1*y1) + dr;
                found_zero = 1;
            }
            x2 = (j-middle_index)*dx;
            y2 = (k-middle_index)*dx;
            u2 = h[j][k];
            j = j + dirs[i].x;
            k = k + dirs[i].y;
        }
        list_R[i] = radius;
        list_tantheta[i] = y1/(x1+dx);  // To avoid division by zero
        //printf("%d\n",i);
    }	
    return 1;
}

int measureRadiusCircularIsland(double**u, int N, double dx, double*x0, double*u0){
	/*	Considers ONLY the section of u(x,y) with y=0
        Finds 4 around one of the two zeros (the one at x>0)
		it returns them in two arrays of 8 elements x0 and u0.
		In python then you can use this data to do a 3rd degree polynomial interpolation
		
	 */
	double central_plateau,x_next,u_next;
    int j = (int)(N/2);
    int i = (int)(N/2);
    central_plateau = u[j][i];
	while (i < N){
	    // Estimate the position of the zero x=xk (y0=0) with a linear fit
		// We store the value of y(=u) until it changes sign, so we can use the previous value to do the fit
		if (u[j][i]*central_plateau < 0){  //as soon as u<0 if the central plateau is >0; as soon as u>0 if the central plateau is <0 
			x0[1] = x_next;
			x0[2] = (i-(int)(N/2))*dx;
			u0[1] = u_next;
			u0[2] = u[j][i];
			//Take other two points around
			x0[0] = ((i-2)-(int)(N/2))*dx;
			u0[0] = u[j][i-2];
			x0[3] = ((i+1)-(int)(N/2))*dx;
			u0[3] = u[j][i+1];
			return 1;
		}
		x_next = (i-(int)(N/2))*dx;
		u_next = u[j][i];
		i = i + 1;
	}
	return 0;
}

double calcellDW(double** hfr, double** hfi, double** q2, int N, double dx){
    double grad2_integral = 0;
    for(int i=0;i<N;i++) {
        for(int j=0;j<N;j++) {
            grad2_integral = grad2_integral + q2[i][j]*(hfr[i][j]*hfr[i][j] + hfi[i][j]*hfi[i][j])*(2*pi/(N*dx))*(2*pi/(N*dx));
        }
	}
	return ((N*dx)*(N*dx))/grad2_integral;
}

double calcDW(double** ghx, double** ghy, double dx, int N){
        /*Compute the L^2/integral of |Grad u|^2. This value, times the thickness of the interface, gives the characteristic lenght \ell_DW*/
        //Remember that FFTW, when performing both DFT and IDFT, does not normalize the data by any factor (nor sqrt(N) nor N in just one of the two)
        double integ_grad2 = 0;
        double L = (double)N*dx;
        for(int i=0;i<N;i++) {
            for(int j=0;j<N;j++) {
                /*We DO NOT take a sqrt for better resolution of the peak position*/
                integ_grad2 = integ_grad2 + (ghx[i][j]*ghx[i][j] + ghy[i][j]*ghy[i][j])*(dx*dx);
            }
        }
        //We do not save the integral of the grad2, but the fraction: TotalArea/integGrad2
        integ_grad2 = (L*L)/integ_grad2;

        return integ_grad2;
}

int calcstructure_fact(double** ufr, double** ufi, int N, double* structure_fac){
    /* Consider S(qx,qy) along the horizontal directions passing through (0, 0) and take the average of the two functions */
    int i,j;
    // Horizontal
    i = 0;
    for (j = 0; j < N; j++){
        structure_fac[j] = ufr[i][j]*ufr[i][j]+ufi[i][j]*ufi[i][j];
    }
    // Vertical
    j = 0;
    for (i = 0; i < N; i++){
        structure_fac[i] = (structure_fac[i] + ufr[i][j]*ufr[i][j]+ufi[i][j]*ufi[i][j])/2;
    }

    return 1;
}

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
        fprintf(dest_file, "%.5f %.20f\n", x, y);
    }
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
		fprintf(dest_file, "%.20f", x);
		for(int j=0; j < leny; j++)
        	fprintf(dest_file, " %.20f", y[j]);
		fprintf(dest_file, "\n");
    }
    fclose(dest_file);
    return 1;
}

int count_lines(FILE* file)
{
    char buf[BUF_SIZE];
    int counter = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    return counter;
}
