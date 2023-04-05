inline void pmf_3d_mult_cells_2(const int *p,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            DG_FP *in_x, DG_FP *in_y,
                            DG_FP *in_z) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int n = 0; n < dg_np; n++) {
    const DG_FP x = in_x[n];
    const DG_FP y = in_y[n];
    const DG_FP z = in_z[n];
    in_x[n] = rx[0] * x + ry[0] * y + rz[0] * z;
    in_y[n] = sx[0] * x + sy[0] * y + sz[0] * z;
    in_z[n] = tx[0] * x + ty[0] * y + tz[0] * z;
  }
}
