inline void poisson_mat_free_mult_cells(const int *p, const double *dr, const double *ds, 
                                        const double *dt, const double *mass, const double *rx,
                                        const double *sx, const double *tx,
                                        const double *ry, const double *sy,
                                        const double *ty, const double *rz,
                                        const double *sz, const double *tz,
                                        const double *J, const double *mm_factor,
                                        const double *in, double *out) {
  const double *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const double *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const double *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const double *mass_mat = &mass[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  
  double tmpX[DG_NP], tmpY[DG_NP], tmpZ[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    out[i] = 0.0;
  }

  for(int m = 0; m < dg_np; m++) {
    double tmpR = 0.0;
    double tmpS = 0.0;
    double tmpT = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      tmpR += dr_mat[ind] * in[n];
      tmpS += ds_mat[ind] * in[n];
      tmpT += dt_mat[ind] * in[n];
    }
    tmpX[m] = rx[0] * tmpR + sx[0] * tmpS + tx[0] * tmpT;
    tmpY[m] = ry[0] * tmpR + sy[0] * tmpS + ty[0] * tmpT;
    tmpZ[m] = rz[0] * tmpR + sz[0] * tmpS + tz[0] * tmpT;
  }

  for(int m = 0; m < dg_np; m++) {
    double tmp0 = 0.0;
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m + n * dg_np;
      tmp0 += mass_mat[ind] * tmpX[n];
      tmp1 += mass_mat[ind] * tmpY[n];
      tmp2 += mass_mat[ind] * tmpZ[n];
    }
    tmpX[m] = J[0] * tmp0;
    tmpY[m] = J[0] * tmp1;
    tmpZ[m] = J[0] * tmp2;
  }

  for(int m = 0; m < dg_np; m++) {
    out[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      int ind = m * dg_np + n;
      out[m] += rx[0] * dr_mat[ind] * tmpX[n];
      out[m] += ry[0] * dr_mat[ind] * tmpY[n];
      out[m] += rz[0] * dr_mat[ind] * tmpZ[n];
      out[m] += sx[0] * ds_mat[ind] * tmpX[n];
      out[m] += sy[0] * ds_mat[ind] * tmpY[n];
      out[m] += sz[0] * ds_mat[ind] * tmpZ[n];
      out[m] += tx[0] * dt_mat[ind] * tmpX[n];
      out[m] += ty[0] * dt_mat[ind] * tmpY[n];
      out[m] += tz[0] * dt_mat[ind] * tmpZ[n];

      out[m] += *mm_factor * J[0] * mass_mat[ind] * in[n];
    }
  }
}