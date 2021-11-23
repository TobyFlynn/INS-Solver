inline void poisson_cubature_op(const double *rx, const double *sx,
                                const double *ry, const double *sy,
                                const double *J, double *op) {
  // Everything in col-major
  double Dx[DG_CUB_NP * DG_NP], Dy[DG_CUB_NP * DG_NP];
  for(int m = 0; m < DG_CUB_NP; m++) {
    // J = -xs.*yr + xr.*ys
    double J_m  = -sx[m] * ry[m] + rx[m] * sy[m];
    // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
    double rx_m = sy[m] / J_m;
    double sx_m = -ry[m] / J_m;
    double ry_m = -sx[m] / J_m;
    double sy_m = rx[m] / J_m;
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_CUB_NP;
      Dx[ind] = rx_m * cubVDr_g[ind] + sx_m * cubVDs_g[ind];
      Dy[ind] = ry_m * cubVDr_g[ind] + sy_m * cubVDs_g[ind];
    }
  }

  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      // op col-major
      int c_ind = m + n * DG_NP;
      // op row-major
      // int c_ind = m * DG_NP + n;
      op[c_ind] = 0.0;
      for(int k = 0; k < DG_CUB_NP; k++) {
        // Dx' and Dy'
        int a_ind = m * DG_CUB_NP + k;
        // Dx and Dy
        int b_ind = n * DG_CUB_NP + k;

        op[c_ind] += J[k] * cubW_g[k] * (Dx[a_ind] * Dx[b_ind] + Dy[a_ind] * Dy[b_ind]);
      }
    }
  }
}
