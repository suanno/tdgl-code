# Initialization

All files in this directory create a folder (if it does not exist) where the simulation will be saved. Inside this folder it creates the files
- state.dat
- init.dat


At the beginning, the files are the same and contain the initial state of the system $u(x, y)$ (first column x, second y and third u).
Header line contains $N$ t $dx$ where the default value of $dx=0.1$ can be changed in the code. 

- datainit: Generates a state with N **random** points uniformly distributed in the range $u_0\pm \epsilon$

        ./datainit N u0 eps foldername

- flatinit: Generates a **flat** state with N points and $u(x)=u_0 \forall x$

        ./flatinit N u0 foldername

- circleinit: Generates a circular domain of radius R and with domain amplitude $\sqrt{C}$: $$u(r)=\sqrt{C}\tanh((r-R)\sqrt{C/2})$$ 

        ./circleinit N C R foldername

- stripeinit: Generates a state with a vertical stripe whose lenght is 1/3 of the system size. The width of interfaces is set by $C$

        ./stripeinit N C foldername

- coffee_grain: Generate a coffe grain subtracting an horizontal stripe of width $d$ from a circular domain of radius R

        ./coffee_grain N C R d foldername