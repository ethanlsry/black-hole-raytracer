#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "boyer-lindquist-metric.h"
#include "ode-solver-rk45-dormand-prince.h"

//test case: unstable spherical photon orbits in Kerr spacetime
void spherical_photon_orbits(){
    double a = 1.0;
    double M = 1.0;

    //activate boyer_lindquist_metric header
    boyer_lindquist_metric metric(a, M);
    
    auto dydx = [&metric] (double x, const std::vector<double> &y) {
        //extract initial conditions
        double r = y[0];
        double theta = y[1];
        double phi = y[2];
        double u_r = y[3];
        double u_theta = y[4];
        double u_phi = y[5];

        //initialize all instace vars
        metric.compute_metric(r, theta);
        double ut_innerRoot = (u_r * u_r * metric.gamma11) + (u_theta * u_theta * metric.gamma22) + (u_phi * u_phi * metric.gamma33);
        double ut = std::sqrt(ut_innerRoot) / metric.alpha; 

        // compute the right hand sides of equations
        std::vector<double> result(6);
        // dr/dt
        result[0] = metric.gamma11 * u_r / ut;

        // dtheta/dt
        result[1] = metric.gamma22 * u_theta / ut;

        // dphi/dt
        result[2] = (metric.gamma33 * u_phi / ut) - metric.beta3;

        // dur/dt
        double res3_term1 = -1 * metric.alpha * ut * metric.d_alpha_dr;
        double res3_term2 = u_phi * metric.d_beta3_dr;
        double res3_term3 = (u_r*u_r*metric.d_gamma11_dr + u_theta*u_theta*metric.d_gamma22_dr + u_phi*u_phi*metric.d_gamma33_dr) / (-2 * ut);
        result[3] = res3_term1 + res3_term2 + res3_term3;

        // dutheta / dt
        double res4_term1 = -1.0 * metric.alpha * ut * metric.d_alpha_dth;
        double res4_term2 = u_phi * metric.d_beta3_dth;
        double res4_term3 = ((u_r*u_r*metric.d_gamma11_dth) + (u_theta*u_theta*metric.d_gamma22_dth) + u_phi*u_phi*metric.d_gamma33_dth) / (-2.0 * ut);
        result[4] = res4_term1 + res4_term2 + res4_term3;

        // duphi / dt
        result[5] = 0;
        
        return result;
    };

    //define initial conditions for spherical photon orbit in Kerr spacetime
    //orbit A
    // double r_initial = 1.0 + std::sqrt(2.0);
    // double theta_initial = M_PI*0.5;
    // double phi_initial = 0;
    // double u_r_initial = 0;
    // double u_theta_initial = std::sqrt(11.0 + (8.0*std::sqrt(2.0)));
    // double u_phi_initial = 0;

    //orbit B
    // double r_initial = 1.0 + std::sqrt(3.0);
    // double theta_initial = M_PI*0.5;
    // double phi_initial = 0;
    // double u_r_initial = 0;
    // double u_theta_initial = std::sqrt(12.0 + 8.0*std::sqrt(3.0));
    // double u_phi_initial = -1.0;

    //orbit C
    // double r_initial = 3.0;
    // double theta_initial = M_PI*0.5;
    // double phi_initial = 0;
    // double u_r_initial = 0;
    // double u_theta_initial = std::sqrt(27.0);
    // double u_phi_initial = -2.0;

    //orbit D
    // double r_initial = 1.0 + 2*std::sqrt(2);
    // double theta_initial = M_PI*0.5;
    // double phi_initial = 0;
    // double u_r_initial = 0;
    // double u_theta_initial = std::sqrt(-13 + 16*std::sqrt(2));
    // double u_phi_initial = -6;
    
    //orbit E
    double r_initial = 2.0;
    double theta_initial = M_PI*0.5;
    double phi_initial = 0;
    double u_r_initial = 0;
    double u_theta_initial = std::sqrt(16);
    double u_phi_initial = 1;

    auto stop = [r_initial](double x, const std::vector<double> &y) { return (y[0]) > (5*r_initial) || (y[0] < (1.2)) || (x>10000); }; 

    metric.compute_metric(r_initial, theta_initial);
    
    std::vector<double> y0 = {
        r_initial, 
        theta_initial, 
        phi_initial,
        u_r_initial,
        u_theta_initial,
        u_phi_initial
    };
    
    double x0 = 0.0;
    double h = 1e-2;

    int num_equations = 6;
    double tolerance_abs = 1e-17;
    double tolerance_rel = 1e-17;

    rk45_dormand_prince rk45(num_equations, tolerance_abs, tolerance_rel);
    rk45.integrate(dydx, stop, h, x0, y0);

    //toggle between A - E depending on desired orbit
    std::ofstream output_spherical("csv-outputs/spherical-photon-output-E.csv");
    for (int i = 0; i < rk45.xs.size(); i++) {
        double r = rk45.result[i][0];
        double theta = rk45.result[i][1];
        double phi = rk45.result[i][2];

        //Pseudo-cartisean coordinates conversion
        double x = std::sqrt(r*r + a*a) * std::sin(theta)*std::cos(phi);
        double y = std::sqrt(r*r + a*a) * std::sin(theta)*std::sin(phi);
        double z = r * std::cos(theta);

        double frac_error = std::abs(r - r_initial) / r_initial;

        output_spherical << rk45.xs[i] << "," << x << "," << y << "," << z << "," << frac_error <<std::endl;
    }
    output_spherical.close();
}