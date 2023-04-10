inline void p_multigrid_vec_minus(const int *order, const DG_FP *vec0,
                                  DG_FP *vec1, const DG_FP *dot) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  for(int i = 0; i < dg_np; i++) {
    vec1[i] -= *dot * vec0[i];
  }
}
