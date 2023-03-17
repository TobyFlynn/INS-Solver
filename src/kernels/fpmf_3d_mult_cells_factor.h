inline void fpmf_3d_mult_cells_factor(const DG_FP *factor,
                                      DG_FP *in_x, DG_FP *in_y, DG_FP *in_z) {
  for(int i = 0; i < DG_NP; i++) {
    in_x[i] *= factor[i];
    in_y[i] *= factor[i];
    in_z[i] *= factor[i];
  }
}
