inline void init_gauss_grad(double *rx, double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {
  // J = -xs.*yr + xr.*ys
  double J[DG_G_NP];
  for(int i = 0; i < DG_G_NP; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < DG_G_NP; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind_row = m * DG_NP + n;
      int ind_col = m + n * DG_GF_NP;
      Dx0[ind_row] = rx[m] * gF0Dr_g[ind_col] + sx[m] * gF0Ds_g[ind_col];
      Dy0[ind_row] = ry[m] * gF0Dr_g[ind_col] + sy[m] * gF0Ds_g[ind_col];
    }
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind_row = m * DG_NP + n;
      int ind_col = m + n * DG_GF_NP;
      Dx1[ind_row] = rx[m + DG_GF_NP] * gF1Dr_g[ind_col] + sx[m + DG_GF_NP] * gF1Ds_g[ind_col];
      Dy1[ind_row] = ry[m + DG_GF_NP] * gF1Dr_g[ind_col] + sy[m + DG_GF_NP] * gF1Ds_g[ind_col];
    }
  }

  for(int m = 0; m < DG_GF_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind_row = m * DG_NP + n;
      int ind_col = m + n * DG_GF_NP;
      Dx2[ind_row] = rx[m + 2 * DG_GF_NP] * gF2Dr_g[ind_col] + sx[m + 2 * DG_GF_NP] * gF2Ds_g[ind_col];
      Dy2[ind_row] = ry[m + 2 * DG_GF_NP] * gF2Dr_g[ind_col] + sy[m + 2 * DG_GF_NP] * gF2Ds_g[ind_col];
    }
  }
}
