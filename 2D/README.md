# Initialization
Compile and run one of the codes in /initialization/ to generate the simulation's folder and the initial state.
See **/initialization/README.md**

Example (random initial conditions):
- Compile

            gcc initialization/datainit.c -fopenmp -lfftw3 -lm -lfftw3_omp -O2 -o .bin/datainit
- Generate the simulation folder and the initial state

              ./datainit N u0 eps foldername

# Generate the C(t) input file
To specify the value of $C(t)$ in
$$\partial_t u = \Delta u + C(t)u -u^3$$
where $u(x,y)$ is the order parameter, you need to generate **fileCin.dat** from **/fileCin generation/generate_fileCin.ipynb**.
It is important to generate this file in the simulation folder (the one specified in the initialization)

# Run simulation
- Compile the simulation code
        
        .gcc tdgl.c observables.c read_write.c -fopenmp -lfftw3 -lm -lfftw3_omp -O2 -o .bin/tdgl

- Run the simulation

        .bin/tdgl tspan foldername

If you run again the last command, the simulation will continue (**not** restart from t=0).

# Observables
Inside the time loop in tdgl.c, you can call the functions defined in **observables.h** to measure different observables as a function of time. 
After the time loop, you can the function **save_observable** to save the observales in a text file.

Remember to allocate (and free) memory for the variables to store this information.