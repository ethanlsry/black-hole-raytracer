#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "boyer-lindquist-metric.h"
#include "ode-solver-rk45-dormand-prince.h"
#include "render-black-hole.h"
//this file renders the black hole image

void render_black_hole(){
    //put ring between r_in = 5M and r_out = 20M
    const double a = 0.99;
    const double M = 1.0;
    const double D = 500.0;
    const double theta_0 = 89.0 * (M_PI / 180.0);

    // const int screenX_dim = 64;
    // const int screenY_dim = 36;
    const int screenX_dim = 1920;
    const int screenY_dim = 1080;

    const int minMax_screenX = 25;
    const int minMax_screenY = 15;

    auto screen = new double[screenY_dim][screenX_dim];

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

    auto stop = [a, M](double x, const std::vector<double> &y) { 
         if (y[0] < 1.01 * (M + std::sqrt(M*M - a*a))) {
            //std::cout << "=";
            return true;
         }
         else if (y[0] > 5*M && y[0] < 20*M && std::abs(y[1] - (0.5 * M_PI))  < 1e-3){
            //std::cout << "*";
            return true;
         }
         else if (y[0] > 1e4){
            //std::cout << "-";
            return true;
         }

         return false;
     };

    //iterate through pixels and calculate trajectory
    for (int y=0; y < screenY_dim; y++){
            double y_sc = ((2.0 * minMax_screenY  * y) / (double) screenY_dim) - minMax_screenY;

            double alpha = y_sc / D;
            double theta_initial = theta_0 - alpha;

        for (int x=0; x < screenX_dim; x++){
            double x_sc = ((2.0 * minMax_screenX  * x) / (double) screenX_dim) - minMax_screenX; 

            double beta = x_sc / D;
            double phi_initial = beta;
            double r_initial = std::sqrt(D*D + x_sc*x_sc + y_sc*y_sc);

            metric.compute_metric(r_initial, theta_initial);
            
            double u_r_initial = -1.0 * std::cos(beta) * std::cos(alpha) * std::sqrt(metric.g_11); 
            double u_theta_initial = std::sin(alpha) * std::sqrt(metric.g_22); 
            double u_phi_initial = std::sin(beta) * std::cos(alpha) * std::sqrt(metric.g_33); 

            std::vector<double> y0 = {
                r_initial, 
                theta_initial, 
                phi_initial,
                u_r_initial,
                u_theta_initial,
                u_phi_initial
            };

            double x0 = 0.0;
            double h = 0.005;

            int num_equations = 6;
            double tolerance_abs = 1e-17; 
            double tolerance_rel = 1e-17; 
            rk45_dormand_prince rk45(num_equations, tolerance_abs, tolerance_rel);

            rk45.integrate(dydx, stop, h, x0, y0);  

            screen[y][x] = 0;

            int lastStep = rk45.xs.size() - 1; 
            double last_r = rk45.result[lastStep][0];
            double last_theta = rk45.result[lastStep][1];

            double last_u_r = 1.0 * rk45.result[lastStep][3];
            double last_u_theta = 1.0 * rk45.result[lastStep][4];
            double last_u_phi = 1.0 * rk45.result[lastStep][5];

            metric.compute_metric(last_r, last_theta);

            double ut_innerRoot = (last_u_r * last_u_r * metric.gamma11) + (last_u_theta * last_u_theta * metric.gamma22) + (last_u_phi * last_u_phi * metric.gamma33);
            double ut = std::sqrt(ut_innerRoot) / metric.alpha; 

            double u_beta_sum = last_u_phi * metric.beta3; // beta1 and beta2 = 0 

            double u_t = -1.0 * metric.alpha * metric.alpha * ut + u_beta_sum; // (eq 16)

            double omega = 1.0 / (a + std::pow(last_r, 1.5)); // M is 1 and omitted
            
            double redshift_factor = (1.0 + (omega * (last_u_phi / ut))) / std::sqrt((-1.0 * metric.g_00) - (omega * omega * metric.g_33) - (2.0 * omega * metric.g_03));
            double intensity = 1.0 / (redshift_factor * redshift_factor * redshift_factor);
            
            if ( last_r > 3*M && last_r < 22*M ){
                screen[y][x] = intensity;
            }
        } 
        std::cout << "(" << y << " of " << screenY_dim << ")" << std::endl;
    }

    std::ofstream output_hole("black-hole-output-HD-theta=89-a=0.99.csv");
    for (int y=0; y < screenY_dim; y++){
        for (int x=0; x < screenX_dim; x++){
            output_hole << std::setprecision(10) << screen[y][x];

            if (x+1 < screenX_dim) {
                output_hole << ",";
            }
        }
        output_hole << std::endl;
    }
    output_hole.close();

    delete[] screen;
}