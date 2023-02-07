inline void fact_poisson_cub_op1(const int *p, const double *cubVDr,
                                 const double *cubVDs, const double *rx,
                                 const double *sx, const double *ry,
                                 const double *sy, const double *J,
                                 const double *factor, double *op) {
  // Get constants
  const int dg_np        = DG_CONSTANTS[(*p - 1) * 5];
  const int dg_cub_np    = DG_CONSTANTS[(*p - 1) * 5 + 2];
  const double *cubVDr_l = &cubVDr[(*p - 1) * DG_CUB_NP * DG_NP];
  const double *cubVDs_l = &cubVDs[(*p - 1) * DG_CUB_NP * DG_NP];
  const double *cubW     = &cubW_g[(*p - 1) * DG_CUB_NP];

  // Everything in col-major
  double Dx[DG_CUB_NP * DG_NP], Dy[DG_CUB_NP * DG_NP];
  for(int m = 0; m < dg_cub_np; m++) {
    // J = -xs.*yr + xr.*ys
    double J_m  = -sx[m] * ry[m] + rx[m] * sy[m];
    // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
    double rx_m = sy[m] / J_m;
    double sx_m = -ry[m] / J_m;
    double ry_m = -sx[m] / J_m;
    double sy_m = rx[m] / J_m;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_cub_np;
      Dx[ind] = rx_m * cubVDr_l[ind] + sx_m * cubVDs_l[ind];
      Dy[ind] = ry_m * cubVDr_l[ind] + sy_m * cubVDs_l[ind];
    }
  }

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      // op col-major
      int c_ind = m + n * dg_np;
      // op row-major
      // int c_ind = m * dg_np + n;
      op[c_ind] = 0.0;
      for(int k = 0; k < dg_cub_np; k++) {
        // Dx' and Dy'
        int a_ind = m * dg_cub_np + k;
        // Dx and Dy
        int b_ind = n * dg_cub_np + k;

        op[c_ind] += J[k] * cubW[k] * factor[k] * (Dx[a_ind] * Dx[b_ind] + Dy[a_ind] * Dy[b_ind]);
        // op[c_ind] += J[k] * cubW[k] * (Dx[a_ind] * Dx[b_ind] + Dy[a_ind] * Dy[b_ind]);
      }
    }
  }
}
