# Cranck-Nicolson scheme in Fourier space for solving the TDGL equation

This code solves
$$\partial_t u = \Delta u + C(t)u -u^3$$
where $u(\mathbf{x})$ is the order parameter and $C(t)$ an input function that depends on time.

For **usage** look in the 1D and 2D directories.

# Required libraries
- FFTW (To compute fourier transforms)
- openmp (For parallel computing)

Parallel computing can be disabled commenting the lines containg the keyword "#pragma".

# Numerical scheme
Crank-Nicolson in Fourier space 
To integrate the TDGL equation, we apply a Fourier transfrom in x, so
$\partial_t u = \partial_{xx}u +C(t)u-\mathcal{F}[u^3]$
becomes ($u(x,t)\rightarrow \mathcal{F}[u(x,t)]=U_q(t)$)
$\partial_t U_q = -q^2U_q+C(t)U_q-\mathcal{F}[\mathcal{F}[u^3]]_q$
So you get rid of the space derivatives and you use the Crank-Nicolson scheme to integrate the equation in time for a small timestep $dt$. Then you do the inverse fourier transform and you retrieve $u(x,t+dt)$. Then you repeat.
## Crank-Nicolson scheme
It is formulated by taking an average of the formulas of Implicit and Explicit schemes:
- Explicit Euler: 

$U(t+dt) = U(t) + [C(t)U(t) - \mathcal{F}[u^3](t)]dt$
- Implicit Euler: 

$U(t+dt) = U(t)+ [C(t+dt)U(t+dt) - \mathcal{F}[u^3](t+dt)]dt$
If we average the two expression (sum them and divide by 2):

$U(t+dt)=U(t)+\frac{dt}2[C(t)U(t)-\mathcal{F}[u^3](t)+C(t+dt)U(t+dt)-\mathcal{F}[u^3](t+dt)]$

Now we make an **approximation** in order to get an explicit formula for $U(t+dt)$ if you know $U(t)$:

$\mathcal{F}[u^3](t+dt)\rightarrow \mathcal{F}[u^3](t)$

After this approximation, we isolate $U(t+dt)$ and we find

$U(t+dt) = U(t)\frac{(1+\frac{dt}{2}C(t))}{(1-\frac{dt}{2}C(t+dt))}-\frac{\mathcal{F}[u^3](t)dt}{(1-\frac{dt}{2}C(t+dt))}$

As you need to compute the FFT of $u^3(t)$ at each step, you need, after each step $dt$, compute the IFFT to get $u(x,t)$, then compute $u^3(x,t)$ and then its FFT. Then you can proceed with the next step.