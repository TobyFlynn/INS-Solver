inline void p_multigrid_vec_norm(const int *order, const double *vec,
                                 double *norm) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    *norm += vec[i] * vec[i];
  }
}
