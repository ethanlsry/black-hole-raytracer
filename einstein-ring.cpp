#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "boyer-lindquist-metric.h"
#include "ode-solver-rk45-dormand-prince.h"
#include "einstein-ring.h"

//this file creates a metric instace and solves b_l system of equations
void einstein_ring(){
    double a = 0.0;
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


    // Case 1 -- Out and back
    int einstein_ring_case = 1;
    auto stop = [](double x, const std::vector<double> &y) { return (y[2] > (2.0 * M_PI) ) ; };
    double r_initial = 10.0 * M;
    double xi = 0.48857874;

    // Case 2 -- Goes around once
    // int einstein_ring_case = 2;
    // auto stop = [](double x, const std::vector<double> &y) { return (y[2] > (4.0 * M_PI) ) ; };
    // double r_initial = 20.0 * M;
    // double xi = 0.24904964;

    double theta_initial = M_PI*0.5;
    double phi_initial = 0;

    metric.compute_metric(r_initial, theta_initial);
    
    double u_r_initial = -1 * std::cos(xi) * std::sqrt(metric.g_11);
    double u_theta_initial = 0;
    double u_phi_initial = std::sin(xi) * std::sqrt(metric.g_33);

    std::vector<double> y0 = {
        r_initial, 
        theta_initial,
        phi_initial,
        u_r_initial,
        u_theta_initial,
        u_phi_initial
    };
    
    double x0 = 0.0;
    double h = 0.01;

    int num_equations = 6;
    double tolerance_abs = 1e-17; 
    double tolerance_rel = 1e-17; 
    rk45_dormand_prince rk45(num_equations, tolerance_abs, tolerance_rel);

    rk45.integrate(dydx, stop, h, x0, y0);

    // Output the x and y coordinates of our orbit
    std::ofstream output_ein("csv-outputs/einstein-ring-output-2.csv");
    for (int i = 0; i < rk45.xs.size(); i++) {
        double r = rk45.result[i][0];
        double theta = rk45.result[i][1];
        double phi = rk45.result[i][2];

        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);

        output_ein << std::setprecision(10) << rk45.xs[i] << "," << x << "," << y << std::endl;
    }
    output_ein.close();


    // Interpolate our last points to find r at phi = 2pi or 4pi
    int lastSpot = rk45.xs.size() - 1;
    
    std::cout << "Points to Interpolate (r,th,phi):\n" 
    << rk45.result[lastSpot-1][0] << ", "
    << rk45.result[lastSpot-1][1] << ", "
    << rk45.result[lastSpot-1][2] << "\n"
    << rk45.result[lastSpot][0] << ", "
    << rk45.result[lastSpot][1] << ", "
    << rk45.result[lastSpot][2] << std::endl;

    double ending_phi = einstein_ring_case == 1 ?  2.0 * M_PI : 4.0 * M_PI;

    std::cout << "\nFind value at phi = " 
    << (einstein_ring_case == 1 ? "2pi" : "4pi")
    << " = " << ending_phi
    << std::endl;

    double delta_r = rk45.result[lastSpot][0] - rk45.result[lastSpot-1][0];
    double delta_phi = rk45.result[lastSpot][2] - rk45.result[lastSpot-1][2];
    
    double slope = delta_r / delta_phi;
    
    double r_estimate = rk45.result[lastSpot-1][0] + ((ending_phi - rk45.result[lastSpot-1][2]) * slope);
    double r_estimate_error = r_estimate - r_initial;
    double r_estimate_error_percent = r_estimate_error / r_initial;

    std::cout << "\nInterpolated r is " << r_estimate 
    << "\nWe expect this value to be " << r_initial
    << "\nThe error is " << r_estimate_error << " or "
    << r_estimate_error_percent << "%"
    << std::endl;
}