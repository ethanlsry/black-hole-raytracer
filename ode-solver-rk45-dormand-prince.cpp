#include <vector>
#include <cmath>
#include <iostream>
#include "ode-solver-rk45-dormand-prince.h"

rk45_dormand_prince::rk45_dormand_prince(int num_equations, double tolerance_abs, double tolerance_rel) {
    // constructor, initialize the data members
    n_eq = num_equations;
    atol = tolerance_abs;
    rtol = tolerance_rel;
    y_tmp.resize(n_eq);
    y_err.resize(n_eq);
    dydx.resize(n_eq);
    dydx_new.resize(n_eq);

    k2.resize(n_eq);
    k3.resize(n_eq);
    k4.resize(n_eq);
    k5.resize(n_eq);
    k6.resize(n_eq);
  }

  
double rk45_dormand_prince::error(const std::vector<double> &y) {
    //compute a scalar error from the error estimate y_err and return it.
    double err = 0.0;
    for (int i = 0; i < n_eq; i++) {
      double scale =
          atol + std::max(std::abs(y[i]), std::abs(y[i])) * rtol;
      err += std::pow(y_err[i] / scale, 2);
      if (std::isnan(err)) {
        std::cout << "err = NaN at i = " << i << "y_err[i]=" << y_err[i] << std::endl;
        exit(1);
      }
    }
    return std::sqrt(err / n_eq);
  };

