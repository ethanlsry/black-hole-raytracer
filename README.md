# Raytracing images of a black hole using General Relativity, C++, and Python
Created in collaboration with Michael Filippini and with support from Alex Chen. Underlying GR theory is from  MacDonald & Thorne (1982), Komissarov (2004), 

This project uses numerical ODE integration to solve the geodesic equation for photons near the event horizon of a rotating black hole. We applied raytracing to render an image of the accretion disk distorted by a black hole's spacetime.

## GR background
A black hole can be described in Kerr spacetime using metric functions which are a mechanism to calcualte the distance between two points separated by d$t$, d$r$, d$\theta$, and d$\phi$. The most basic parameters in these functions are M, mass of the black hole, and a, the spin parameter. In this project we use the Boyer-Lindquist formulation of metric functions.

## Einstein Ring in Schwarzschild Spacetime
Our first test case was plotting Einstein rings in Schwarzchild spacetime. This test was used to determine the error from our RK45 solver as well as to confirm that our metric and differential equations were correct. It is worth noting that because this black hole was not rotating (a = 0), some errors in the metric or differential equation specifications were not caught until our next test case.
The test case was run by calling einstein_ring() in einstein-ring.cpp. This file creates a metric instance, the six partial differential equations which describe trajectories near a black hole, and its stop condition. Initial conditions for both test cases are included, but one should always remain commented. After running the integration, the last two points are used to calculate r at either $\phi$ = 2$\pi$or 4$\pi$.
The first case for Einstein rings was a photon that travels behind the black hole and then back to its starting position. For this case, we will use r = 10M and critical angle ξ = 0.48857874. Critical angle values were given by Müller (2008). We created a metric instance inside einstein ring.cpp. Figure 1 was generated with RK45 absolute and relative tolerances set to e−17.

![FIG. 1 Einstein Ring Case 1](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/einstein-ring-output-case-1.png?raw=true)

The second case for Einstein rings was a photon that travels one time around the black hole and then back to its starting position. For this case, we used r = 20M and ξ = 0.24904964. Figure 2 was generated with RK45 absolute and relative tolerances set to e−17.

![FIG.21 Einstein Ring Case 2](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/einstein-ring-output-case-2.png?raw=true)


## File Structure
Our files are set up such that they can all be compiled together without conflicts. We only have a main method in raytracer-runner.cpp and we properly use header files for declarations. This meant changing the provided Dormand Prince solver to keep only template functions and declarations in a header file and function definitions in a cpp file. Without this change we ran the risk of multiple definitions of functions when including the Dormand Prince solver in multiple files. In future developemnt we hope to render multiple consecutive initial parameters to create a black hole movie.