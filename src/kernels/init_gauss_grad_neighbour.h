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
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx0[ind] = rx[i] * gF0DrR_g[ind] + sx[i] * gF0DsR_g[ind];
        Dy0[ind] = ry[i] * gF0DrR_g[ind] + sy[i] * gF0DsR_g[ind];
      }
    }
  } else {
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx0[ind] = rx[i] * gF0Dr_g[ind] + sx[i] * gF0Ds_g[ind];
        Dy0[ind] = ry[i] * gF0Dr_g[ind] + sy[i] * gF0Ds_g[ind];
      }
    }
  }

  if(reverse[1]) {
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx1[ind] = rx[i + DG_GF_NP] * gF1DrR_g[ind] + sx[i + DG_GF_NP] * gF1DsR_g[ind];
        Dy1[ind] = ry[i + DG_GF_NP] * gF1DrR_g[ind] + sy[i + DG_GF_NP] * gF1DsR_g[ind];
      }
    }
  } else {
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx1[ind] = rx[i + DG_GF_NP] * gF1Dr_g[ind] + sx[i + DG_GF_NP] * gF1Ds_g[ind];
        Dy1[ind] = ry[i + DG_GF_NP] * gF1Dr_g[ind] + sy[i + DG_GF_NP] * gF1Ds_g[ind];
      }
    }
  }

  if(reverse[2]) {
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx2[ind] = rx[i + 2 * DG_GF_NP] * gF2DrR_g[ind] + sx[i + 2 * DG_GF_NP] * gF2DsR_g[ind];
        Dy2[ind] = ry[i + 2 * DG_GF_NP] * gF2DrR_g[ind] + sy[i + 2 * DG_GF_NP] * gF2DsR_g[ind];
      }
    }
  } else {
    for(int j = 0; j < DG_NP; j++) {
      for(int i = 0; i < DG_GF_NP; i++) {
        int ind = j * DG_GF_NP + i;
        Dx2[ind] = rx[i + 2 * DG_GF_NP] * gF2Dr_g[ind] + sx[i + 2 * DG_GF_NP] * gF2Ds_g[ind];
        Dy2[ind] = ry[i + 2 * DG_GF_NP] * gF2Dr_g[ind] + sy[i + 2 * DG_GF_NP] * gF2Ds_g[ind];
      }
    }
  }
}
