inline void euler_l2_vortex_error_0(const double *time, const double *x,
                                    const double *y, const double *J,
                                    double *rho) {
  const double centre_x = 0.0;
  const double centre_y = fmod(*time, 40.0) > 20.0 ? fmod(*time, 40.0) - 40.0 : fmod(*time, 40.0);
  const double PI = 3.141592653589793238463;
  const double R = 1.5;
  const double S = 13.5;
  const double M = 0.4;
  const double *cubW  = &cubW_g[(DG_ORDER - 1) * DG_CUB_NP];

  for(int i = 0; i < DG_CUB_NP; i++) {
    double rel_x = x[i] - centre_x;
    double rel_y = y[i] - centre_y;
    if(centre_y > 18.0 && y[i] < -18.0) {
      rel_y = y[i] - (centre_y - 40.0);
    }
    if(centre_y < -18.0 && y[i] > 18.0) {
      rel_y = y[i] - (centre_y + 40.0);
    }

    double f = (1.0 - rel_x * rel_x - rel_y * rel_y) / (2.0 * R * R);
    double rho_tmp = 1.0 - ((S * S * M * M * (gamma_e - 1.0) * exp(2.0 * f)) / (8.0 * PI * PI));
    double rho_ = pow(rho_tmp, 1.0 / (gamma_e - 1.0));

    if(fabs(rel_x) <= 2.0 && fabs(rel_y) <= 2.0)
      rho[i] = cubW[i] * J[i] * (rho[i] - rho_) * (rho[i] - rho_);
    else
      rho[i] = 0.0;
  }
}
