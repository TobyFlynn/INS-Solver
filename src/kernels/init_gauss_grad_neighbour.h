inline void init_gauss_grad_neighbour(const int *reverse, double *rx,
                            double *sx, double *ry,  double *sy,
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

  if(reverse[0]) {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx0[ind] = rx[m] * gF0DrR_g[ind] + sx[m] * gF0DsR_g[ind];
        Dy0[ind] = ry[m] * gF0DrR_g[ind] + sy[m] * gF0DsR_g[ind];
      }
    }
  } else {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx0[ind] = rx[m] * gF0Dr_g[ind] + sx[m] * gF0Ds_g[ind];
        Dy0[ind] = ry[m] * gF0Dr_g[ind] + sy[m] * gF0Ds_g[ind];
      }
    }
  }

  if(reverse[1]) {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx1[ind] = rx[m + DG_GF_NP] * gF1DrR_g[ind] + sx[m + DG_GF_NP] * gF1DsR_g[ind];
        Dy1[ind] = ry[m + DG_GF_NP] * gF1DrR_g[ind] + sy[m + DG_GF_NP] * gF1DsR_g[ind];
      }
    }
  } else {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx1[ind] = rx[m + DG_GF_NP] * gF1Dr_g[ind] + sx[m + DG_GF_NP] * gF1Ds_g[ind];
        Dy1[ind] = ry[m + DG_GF_NP] * gF1Dr_g[ind] + sy[m + DG_GF_NP] * gF1Ds_g[ind];
      }
    }
  }

  if(reverse[2]) {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx2[ind] = rx[m + 2 * DG_GF_NP] * gF2DrR_g[ind] + sx[m + 2 * DG_GF_NP] * gF2DsR_g[ind];
        Dy2[ind] = ry[m + 2 * DG_GF_NP] * gF2DrR_g[ind] + sy[m + 2 * DG_GF_NP] * gF2DsR_g[ind];
      }
    }
  } else {
    for(int m = 0; m < DG_GF_NP; m++) {
      for(int n = 0; n < DG_NP; n++) {
        int ind = m * DG_NP + n;
        Dx2[ind] = rx[m + 2 * DG_GF_NP] * gF2Dr_g[ind] + sx[m + 2 * DG_GF_NP] * gF2Ds_g[ind];
        Dy2[ind] = ry[m + 2 * DG_GF_NP] * gF2Dr_g[ind] + sy[m + 2 * DG_GF_NP] * gF2Ds_g[ind];
      }
    }
  }
}
