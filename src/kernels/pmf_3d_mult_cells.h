inline void pmf_3d_mult_cells(const int *p, const DG_FP *dr,
                            const DG_FP *ds, const DG_FP *dt, const DG_FP *mass,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *J, const DG_FP *l_x, const DG_FP *l_y,
                            const DG_FP *l_z, const DG_FP *in_x, const DG_FP *in_y,
                            const DG_FP *in_z, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *mass_mat = &mass[(*p - 1) * DG_NP * DG_NP];
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  DG_FP tmpX[DG_NP], tmpY[DG_NP], tmpZ[DG_NP];

  for(int m = 0; m < dg_np; m++) {
    tmpX[m] = 0.0;
    tmpY[m] = 0.0;
    tmpZ[m] = 0.0;
    for(int n = 0; n < dg_np; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, dg_np, dg_np);
      tmpX[m] += mass_mat[ind] * in_x[n];
      tmpY[m] += mass_mat[ind] * in_y[n];
      tmpZ[m] += mass_mat[ind] * in_z[n];
    }
    tmpX[m] *= J[0];
    tmpY[m] *= J[0];
    tmpZ[m] *= J[0];
  }

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < dg_np; n++) {
      // int ind = m * dg_np + n;
      int ind = DG_MAT_IND(n, m, dg_np, dg_np);
      out[m] += dr_mat[ind] * (rx[0] * (tmpX[n] + l_x[n]) + ry[0] * (tmpY[n] + l_y[n]) + rz[0] * (tmpZ[n] + l_z[n]));
      out[m] += ds_mat[ind] * (sx[0] * (tmpX[n] + l_x[n]) + sy[0] * (tmpY[n] + l_y[n]) + sz[0] * (tmpZ[n] + l_z[n]));
      out[m] += dt_mat[ind] * (tx[0] * (tmpX[n] + l_x[n]) + ty[0] * (tmpY[n] + l_y[n]) + tz[0] * (tmpZ[n] + l_z[n]));
    }
  }
}
