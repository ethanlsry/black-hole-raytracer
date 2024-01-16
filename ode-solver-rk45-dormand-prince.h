#include <vector>

//credit to Alex Chen
class rk45_dormand_prince {
public:
    rk45_dormand_prince(int num_equations, double tolerance_abs, double tolerance_rel);

    template <typename F, typename StopCondition>
    std::vector<double> integrate(const F &f, const StopCondition &stop, double h,
                                double x0, const std::vector<double> &y0) {
    // clear the arrays so that we start fresh every time we call this function
    xs.clear();
    result.clear();

    // initial conditions
    double x = x0;
    std::vector<double> y = y0;
    double err = 1.0;
    double err_prev = 1.0;
    dydx = f(x, y);

    // Always push the initial condition into the results
    xs.push_back(x);
    result.push_back(y);

    while (stop(x, y) == false) {
      y = step(f, h, x, y);
      err = error(y);
      // If err is fortuitously too small, set it to some lower bound
      err = std::max(err, 1.0e-10);

      // Accept the step if the scalar error is below 1, otherwise reject it and
      // do not move forward
      if (err < 1.0) {
        x += h;
        xs.push_back(x);
        result.push_back(y);
        dydx = dydx_new;
      }

      // Adjust h as needed
      double err_alpha = 0.7 / 5.0;
      double err_beta = 0.4 / 5.0;
      h = std::max(hmin, 0.9 * h * std::pow(err, -err_alpha) *
                             std::pow(err_prev, err_beta));
      h = std::min(hmax, h);
      err_prev = err;
    }
    return y;
  };

    
  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    // Compute the next step in y, given x and y of the current step
    std::vector<double> y_next(n_eq);

    // First step
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * a21 * dydx[i];
    }

    // Second step
    k2 = f(x + c2 * h, y_tmp);
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
    }

    // Third step
    k3 = f(x + c3 * h, y_tmp);
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i]);
    }

    // Fourth step
    k4 = f(x + c4 * h, y_tmp);
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i] +
                              a54 * k4[i]);
    }

    // Fifth step
    k5 = f(x + c5 * h, y_tmp);
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i] +
                              a64 * k4[i] + a65 * k5[i]);
    }

    // Sixth step
    k6 = f(x + h, y_tmp);
    for (int i = 0; i < n_eq; i++) {
      y_next[i] = y[i] + h * (a71 * dydx[i] + a72 * k2[i] + a73 * k3[i] +
                              a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    }
    dydx_new = f(x + h, y_next);

    // Estimate y_err for each y in the vector using the difference
    // between y1 and y2
    for (int i = 0; i < n_eq; i++) {
      y_err[i] = h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] +
                       e6 * k6[i] + e7 * dydx_new[i]);
    }

    return y_next;
  };

    
  double error(const std::vector<double> &y);

    int n_eq;
    double atol, rtol;
    double hmin = 1.0e-15;
    double hmax = 1.0;
    std::vector<double> dydx, dydx_new;
    std::vector<double> k2, k3, k4, k5, k6, y_tmp, y_err;
    std::vector<double> xs;
    std::vector<std::vector<double>> result;

    const double c2 = 1.0 / 5.0;
    const double c3 = 3.0 / 10.0;
    const double c4 = 4.0 / 5.0;
    const double c5 = 8.0 / 9.0;

    const double a21 = 1.0 / 5.0;
    const double a31 = 3.0 / 40.0;
    const double a32 = 9.0 / 40.0;
    const double a41 = 44.0 / 45.0;
    const double a42 = -56.0 / 15.0;
    const double a43 = 32.0 / 9.0;
    const double a51 = 19372.0 / 6561.0;
    const double a52 = -25360.0 / 2187.0;
    const double a53 = 64448.0 / 6561.0;
    const double a54 = -212.0 / 729.0;
    const double a61 = 9017.0 / 3168.0;
    const double a62 = -355.0 / 33.0;
    const double a63 = 46732.0 / 5247.0;
    const double a64 = 49.0 / 176.0;
    const double a65 = -5103.0 / 18656.0;
    const double a71 = 35.0 / 384.0;
    const double a72 = 0.0;
    const double a73 = 500.0 / 1113.0;
    const double a74 = 125.0 / 192.0;
    const double a75 = -2187.0 / 6784.0;
    const double a76 = 11.0 / 84.0;

    const double e1 = 71.0 / 57600.0;
    const double e2 = 0.0;
    const double e3 = -71.0 / 16695.0;
    const double e4 = 71.0 / 1920.0;
    const double e5 = -17253.0 / 339200.0;
    const double e6 = 22.0 / 525.0;
    const double e7 = -1.0 / 40.0;

    const double d1 = -12715105075.0 / 11282082432.0;
    const double d2 = 0.0;
    const double d3 = 87487479700.0 / 32700410799.0;
    const double d4 = -10690763975.0 / 1880347072.0;
    const double d5 = 701980252875.0 / 199316789632.0;
    const double d6 = -1453857185.0 / 822651844.0;
    const double d7 = 69997945.0 / 29380423.0;
};