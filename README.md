# Raytracing images of a black hole using General Relativity, C++, and Python
Created in collaboration with Michael Filippini and with support from Alex Chen. Underlying GR theory is from  MacDonald & Thorne (1982), Komissarov (2004), 

This project uses numerical ODE integration to solve the geodesic equation for photons near the event horizon of a rotating black hole. We applied raytracing to render an image of the accretion disk distorted by a black hole's spacetime.

## GR background
A black hole can be described in Kerr spacetime using metric functions which are a mechanism to calcualte the distance between two points separated by d$t$, d$r$, d$\theta$, and d$\phi$. The most basic parameters in these functions are M, mass of the black hole, and a, the spin parameter. In this project we use the Boyer-Lindquist formulation of metric functions.

## Einstein Ring in Schwarzschild Spacetime
Our first test case was plotting Einstein rings in Schwarzchild spacetime. This test was used to determine the error from our RK45 solver as well as to confirm that our metric and differential equations were correct. It is worth noting that because this black hole was not rotating (a = 0), some errors in the metric or differential equation specifications were not caught until our next test case.
The test case was run by calling einstein_ring() in einstein-ring.cpp. This file creates a metric instance, the six partial differential equations which describe trajectories near a black hole, and its stop condition. Initial conditions for both test cases are included, but one should always remain commented. After running the integration, the last two points are used to calculate r at either $\phi$ = 2$\pi$or 4$\pi$.
The first case for Einstein rings was a photon that travels behind the black hole and then back to its starting position. For this case, we will use r = 10M and critical angle ξ = 0.48857874. Critical angle values were given by Müller (2008). We created a metric instance inside einstein ring.cpp. Figure 1 was generated with RK45 absolute and relative tolerances set to 1.0e−17.

![FIG. 1 Einstein Ring Case 1](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/einstein-ring-output-case-1.png?raw=true)

The second case for Einstein rings was a photon that travels one time around the black hole and then back to its starting position. For this case, we used r = 20M and ξ = 0.24904964. Figure 2 was generated with RK45 absolute and relative tolerances set to 1.0e−17.

![FIG.2 Einstein Ring Case 2](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/einstein-ring-output-case-2.png?raw=true)


## Unstable Spherical Photon Orbits in Kerr Spacetime

We implemented five unstable photon orbits with different initial parameters in spherical-photon-orbits.cpp. 3D plots (via K3D) showed each trajectory continued for some amount of time before it either fell to the singularity or its radius grew large. We decided to define “large” in our stop function as five times the initial radius. Orbits A, C, D, and E eventually spiraled outward while orbit B eventually collapsed inward. Orbit D is shown below (FIG 3).

![FIG.3: 3D Plots of Unstable Photon Orbit](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/spherical-orbit-k3d-plot.png?raw=true)

We also calculated fractional error of radii as functions of time to compare stability across orbits. We found that orbit A was the most stable while orbit B was the least stable based on duration in orbit (FIG 4).

![FIG. 4: Decay of Unstable Orbits](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/spherical_photon_orbit_frac_error.png?raw=true)

These test cases established strong confidence in our boyer lindquist metric file and lambda function. They worked for not only the Einstein rings, but also a scenario in which the a parameter representing rotation of the black hole is set to a non-zero value. Now we can proceed to the fun part: raytracing a black hole.


## Rendering a Black Hole Image
After validating our geodesic equation solutions for photons near a black hole, we built render-black-hole.cpp for raytracing. There are two coordinate systems we use to produce the images. The array holding the brightness values has integer coordinates (indices) that determine the resolution of the image. There is also a coordinate system determining the size of the screen relative to the black hole. This determines the frame of the image and it always set with xscreen ∈ [−25,25] and yscreen ∈ [−15,15] to capture the entire black hole. These coordinate systems are defined in the first section of render black hole.cpp. The translation between these coordinate systems is a simple linear map necessary to determine the initial conditions of a photon whose value goes in a specific index of the array.
The stop function for render-black-hole.cpp has three cases: the photon is pulled into the black hole, the photon hits the accretion disk, or the photon flies far away from the black hole. If the photon’s radius to the black hole is within 1% of the event horizon's radius it falls in to the black hole. If the photon’s radius is within [5M, 20M ] of the black hole and $\theta$ is within 1.0e−3 of π/2, the photon hit the accretion disk. If the radius is greater than 1.0e4, it flies away from the black hole. These stop conditions and tolerances balance the accuracy of the image
and the computational cost to create the image.
Once the integration stops, we use the final radius to determine if it hit the accretion disk. If so, we can then
calculate the total redshift factor to determine its intensity. If not, it will have an intensity of 0 and show up black in our image.
We iterated over each pixel based on a desired resolution and solved the geodesic equation to determine whether or not a photon would reach the observer at this point. In the case of a photon presence calculated the scalar light intensity value based on the effect of gravitational redshift and Doppler boost.[3] [4] We used constant M = 1.0 and observer distance D = 500 while varying rotation variable a and inclination θ0. Initially we rendered low resolution images to gain confidence in our output before committing to multiple hours of rendering full HD images. FIG 5 shows a low definition black hole image using a = 0.99 and $\theta$ = 85 $\deg$ which took around 30 seconds on a M1 processor. FIG 6 shows a HD image using a = 0.99 and $\theta$ = 85 $\deg$ which took around 8 hours on a M1 processor. The black hole shadow, photon ring, and accretion disk are all visible. The left side of the accretion disk is brighter than the right side due to the effects of Doppler beaming from black hole rotation. It is possible to increase the parameters to render black hole images with different spin parameters and inclinations.

![FIG. 5: Ray Traced Black Hole (64 × 36 resolution, a = 0.99, inclination $\theta$ = 85$\deg$)](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/ray-traced-black-hole-64-36.png?raw=true)

![FIG. 6: HD Ray Traced Black Hole (1920 × 1080 resolution, a = 0.99, inclination $\theta$ = 85$\deg$)](https://github.com/ethanlsry/black-hole-raytracer/blob/main/plots/ray-traced-black-hole-1920-1080.png?raw=true)

## File Structure
Our files are set up such that they can all be compiled together without conflicts. We only have a main method in raytracer-runner.cpp and we properly use header files for declarations. This meant changing the provided Dormand Prince solver to keep only template functions and declarations in a header file and function definitions in a cpp file. Without this change we ran the risk of multiple definitions of functions when including the Dormand Prince solver in multiple files. In future developemnt we hope to render multiple consecutive initial parameters to create a black hole movie.