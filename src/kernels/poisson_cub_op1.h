inline void poisson_cub_op1(const int *p, const DG_FP *cubVDr,
                            const DG_FP *cubVDs, const DG_FP *rx,
                            const DG_FP *sx, const DG_FP *ry,
                            const DG_FP *sy, const DG_FP *J,
                            DG_FP *op) {
  // Get constants
  const int dg_np        = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_cub_np    = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 2];
  const DG_FP *cubVDr_l = &cubVDr[(*p - 1) * DG_CUB_NP * DG_NP];
  const DG_FP *cubVDs_l = &cubVDs[(*p - 1) * DG_CUB_NP * DG_NP];
  const DG_FP *cubW     = &cubW_g[(*p - 1) * DG_CUB_NP];

  // Everything in col-major
  DG_FP Dx[DG_CUB_NP * DG_NP], Dy[DG_CUB_NP * DG_NP];
  for(int m = 0; m < dg_cub_np; m++) {
    // J = -xs.*yr + xr.*ys
    DG_FP J_m  = -sx[m] * ry[m] + rx[m] * sy[m];
    // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
    DG_FP rx_m = sy[m] / J_m;
    DG_FP sx_m = -ry[m] / J_m;
    DG_FP ry_m = -sx[m] / J_m;
    DG_FP sy_m = rx[m] / J_m;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m + n * dg_cub_np;
      int ind = DG_MAT_IND(m, n, dg_cub_np, dg_np);
      Dx[ind] = rx_m * cubVDr_l[ind] + sx_m * cubVDs_l[ind];
      Dy[ind] = ry_m * cubVDr_l[ind] + sy_m * cubVDs_l[ind];
    }
  }

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      // op col-major
      // int c_ind = m + n * dg_np;
      int c_ind = DG_MAT_IND(m, n, dg_np, dg_np);
      // op row-major
      // int c_ind = m * dg_np + n;
      op[c_ind] = 0.0;
      for(int k = 0; k < dg_cub_np; k++) {
        // Dx' and Dy'
        // int a_ind = m * dg_cub_np + k;
        int a_ind = DG_MAT_IND(k, m, dg_cub_np, dg_np);
        // Dx and Dy
        // int b_ind = n * dg_cub_np + k;
        int b_ind = DG_MAT_IND(k, n, dg_cub_np, dg_np);

        op[c_ind] += J[k] * cubW[k] * (Dx[a_ind] * Dx[b_ind] + Dy[a_ind] * Dy[b_ind]);
        // op[c_ind] += J[k] * cubW[k] * (Dx[a_ind] * Dx[b_ind] + Dy[a_ind] * Dy[b_ind]);
      }
    }
  }
}
