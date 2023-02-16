inline void poisson_matrix_free_3d_mm_mult_cells(const int *p, const DG_FP *dr, 
                            const DG_FP *ds, const DG_FP *dt, const DG_FP *mass, 
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty, 
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *J, const DG_FP *mm_factor, 
                            const DG_FP *in, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  DG_FP tmpX[DG_NP], tmpY[DG_NP], tmpZ[DG_NP];
  for(int m = 0; m < dg_np; m++) {
    DG_FP tmpR = 0.0;
    DG_FP tmpS = 0.0;
    DG_FP tmpT = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      tmpR += dr_mat[ind] * in[n];
      tmpS += ds_mat[ind] * in[n];
      tmpT += dt_mat[ind] * in[n];
    }
    tmpX[m] = rx[0] * tmpR + sx[0] * tmpS + tx[0] * tmpT;
    tmpY[m] = ry[0] * tmpR + sy[0] * tmpS + ty[0] * tmpT;
    tmpZ[m] = rz[0] * tmpR + sz[0] * tmpS + tz[0] * tmpT;
  }

  for(int m = 0; m < dg_np; m++) {
    DG_FP tmp0 = 0.0;
    DG_FP tmp1 = 0.0;
    DG_FP tmp2 = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
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
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(n, m, dg_np, dg_np);
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
