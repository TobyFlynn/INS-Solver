inline void fpmf_3d_mult_cells_factor(const int *p, const DG_FP *factor,
                                      DG_FP *in_x, DG_FP *in_y, DG_FP *in_z) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    in_x[i] *= factor[i];
    in_y[i] *= factor[i];
    in_z[i] *= factor[i];
  }
}
