class boyer_lindquist_metric {
public:

    //declare all instance vars
    double a, M;
    double alpha, beta3;
    double gamma11, gamma22, gamma33; //upper gamma^(ij) components
    double g_00, g_03, g_11, g_22, g_33; //lower g_(mu nu) components
    double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr;
    double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;

    boyer_lindquist_metric (double a0, double M0);

    void compute_metric (double r, double theta);
};