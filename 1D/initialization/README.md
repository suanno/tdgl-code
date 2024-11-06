# Initialization

**ATTENTION**: Please remove the "else" cases for the terminal inputs. They are dangerous, it's better to receive an error!

All the following codes generate a file "fileinit.dat". This file is not sufficient for running "tdgl.c", but a previous run of "tdglfd.c" is necessary to generate the file "tdgl_result.dat" from the information contained in "fileinit.dat".
If you run "tdglfd.c" with $\Delta t = 0$, then the content of "fileinit.dat" is simply copyed in "tdgl_result.dat".


- datainit: Generates a state with N **random** points randomly distributed around $u_0$ with an amplitude $\pm$ eps

        ./datainit N eps u0

- flatinit: Generates a **flat** state with N points

        ./flatinit N u0

- flateps: Generates a **flat** profile $u_0$ of N point _with added_ a small perturbation.
The perturbation consist in giving a small amplitude eps to the modes k=1,2,...,N/2 (the others have the same period, but they are "moving" in opposite direction, so don't care) and randomly dephasing these modes (to avoid particular destructive/constructive interference)

        ./flateps N u0 eps

- flatsinit: Generates a **flat** profile $u_0$ of N points to which is added a sine wave with wavelenght $\lambda$ and amplitude eps.
    
    In order to enable the choice of $\lambda$ you need to provide $L$ and $dx$.
    For this reason, executing this script will immediately create "tdgl_result.dat", without the need of executing "tdglfd" as usual.
    Here is not necessary to smoother the profile, at it is already so. So there aren't any _numerical_ issues if we do not evolve (for the first steps) the state with explicit euler in real space. 

        ./flatsinit L dx u0 lamb eps

- twokinksinit: Generate a **two kinks profile** with tanh shape

        # To implement...