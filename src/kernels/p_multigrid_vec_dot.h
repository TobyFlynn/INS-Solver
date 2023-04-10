inline void p_multigrid_vec_dot(const int *order, const DG_FP *vec0,
                                const DG_FP *vec1, DG_FP *dot) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS];

  for(int i = 0; i < dg_np; i++) {
    *dot += vec0[i] * vec1[i];
  }
}
