# Initialization

All files in this directory create a folder (if it does not exist) where the simulation will be saved. Inside this folder it creates the files
- state.dat
- init.dat


At the beginning, the files are the same and contain the initial state of the system $u(x)$ (first column is x, second is u).
Header line contains $N$ t $dx$ where the default value of $dx=0.1$ can be changed in the code. 

- datainit: Generates a state with N **random** points uniformly distributed in the range $u_0\pm \epsilon$

        ./datainit N u0 eps foldername

- flatinit: Generates a **flat** state with N points and $u(x)=u_0 \forall x$

        ./flatinit N u0 foldername


- flatsinit: Generates a **flat** profile $u_0$ of N points to which is added a sine wave with wavelenght $\lambda$ and amplitude eps.

        ./flatsinit N dx u0 lamb eps foldername

- twokinksinit: Generate a state with a kink followed by an anti-kink at distance d. The amplitude of the domains is u0 and the shape of the kink/anti-kink is tanh

        ./twokinksinit N u0 d foldername