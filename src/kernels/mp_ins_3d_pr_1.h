inline void mp_ins_3d_pr_1(const double *b0, const double *b1, const double *dt,
                           const double *dPdN, double *dPdN_old, double *divVelT,
                           const double *rho, double *pr_factor) {
  double factor = 1.0 / *dt;
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor;
    pr_factor[i] = 1.0 / rho[i];
  }

  for(int i = 0; i < 4 * DG_NPF; i++) {
    dPdN_old[i] = (*b0) * dPdN[i] + (*b1) * dPdN_old[i];
  }
}