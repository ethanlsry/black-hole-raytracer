#include "boyer-lindquist-metric.h"
#include <cmath>

boyer_lindquist_metric::boyer_lindquist_metric (double a0, double M0){
    //initialize instance vars
    a = a0;
    M = M0;
}

void boyer_lindquist_metric::compute_metric (double r, double theta) {
    //pre-calculate sin and cos
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);

    //compute Kerr spacetime metric functions parameters
    double rho_squared = r*r + a*a*cos_theta*cos_theta;
    double delta = r*r + a*a - 2*M*r;
    double summation_symbol = (r*r + a*a) * (r*r + a*a) - a*a*delta*sin_theta*sin_theta;

    //compute lapse function
    this->alpha = std::sqrt(rho_squared * delta / summation_symbol);

    //compute upper index shift vector
    this->beta3 = (-2*M*a*r / summation_symbol);

    //components of upper gamma
    this->gamma11 = delta / rho_squared;
    this->gamma22 = 1 / rho_squared;
    this->gamma33 = rho_squared / (summation_symbol * sin_theta * sin_theta);

    //components of lower gamma
    this->g_00 = ((2 * M * r) / rho_squared) - 1;
    this->g_03 = ((-2 * M * a * r * sin_theta * sin_theta) / rho_squared);
    this->g_11 = rho_squared / delta;
    this->g_22 = rho_squared;
    this->g_33 = summation_symbol * sin_theta * sin_theta / rho_squared;

    //derivatives wrt r:
    double d_rhoSquared_dr = 2.0 * r;
    double d_delta_dr = 2.0 * r - 2.0 * M;
    double d_summation_dr = (2.0 * (r * r + a * a) * 2.0 * r) - (a * a * sin_theta * sin_theta * d_delta_dr);

    double d_alpha_dr_term1 = rho_squared * delta * ((-1.0 * d_summation_dr)/(summation_symbol * summation_symbol));
    double d_alpha_dr_term2 = (rho_squared / summation_symbol) * d_delta_dr; 
    double d_alpha_dr_term3 = (delta / summation_symbol) * d_rhoSquared_dr;

    this->d_alpha_dr = (0.5 * std::sqrt(summation_symbol / (rho_squared * delta))) * (d_alpha_dr_term1 + d_alpha_dr_term2 + d_alpha_dr_term3);
    this->d_beta3_dr = (-2.0 * M * a) * ((1.0 / summation_symbol) + ((-1.0 * d_summation_dr * r)/(summation_symbol * summation_symbol)));

    this->d_gamma11_dr = (d_delta_dr / rho_squared) + (delta * (-1.0 * d_rhoSquared_dr) / (rho_squared * rho_squared));
    this->d_gamma22_dr = (-1.0 * d_rhoSquared_dr) / (rho_squared * rho_squared);
    this->d_gamma33_dr = (d_rhoSquared_dr / (summation_symbol*sin_theta*sin_theta)) + ((rho_squared/(sin_theta*sin_theta)) * ((-1.0 * d_summation_dr)/(summation_symbol * summation_symbol)));

    //derivatives wrt theta:
    double d_rhoSquared_dth = -2 * a * a * sin_theta * cos_theta;
    double d_sinSquared_dth = 2 * cos_theta * sin_theta; 
    double d_summation_dth = -1 * a * a * delta * d_sinSquared_dth; 

    double d_alpha_dth_term1 = d_rhoSquared_dth * (delta / summation_symbol);
    double d_alpha_dth_term2 = ((-1.0 * d_summation_dth) / (summation_symbol * summation_symbol)) * (rho_squared * delta);
   
    this->d_alpha_dth = (0.5 * std::sqrt(summation_symbol / (rho_squared * delta))) * (d_alpha_dth_term1 + d_alpha_dth_term2);     
    this->d_beta3_dth = (-2.0 * M * a * r) * ((-1.0 * d_summation_dth) / (summation_symbol * summation_symbol));

    this->d_gamma11_dth = (-1.0*delta*d_rhoSquared_dth) / (rho_squared*rho_squared);
    this->d_gamma22_dth =  (-1.0*d_rhoSquared_dth) / (rho_squared*rho_squared); 

    double d_gamma33_dth_term1 = d_rhoSquared_dth / (summation_symbol * sin_theta * sin_theta);
    double d_gamma33_dth_term2 = (rho_squared/summation_symbol) * ((-1 * d_sinSquared_dth ) / (sin_theta * sin_theta * sin_theta * sin_theta));
    double d_gamma33_dth_term3 = (rho_squared/(sin_theta * sin_theta)) * ((-1 * d_summation_dth) / (summation_symbol * summation_symbol));

    this->d_gamma33_dth = d_gamma33_dth_term1 + d_gamma33_dth_term2 + d_gamma33_dth_term3; 
}