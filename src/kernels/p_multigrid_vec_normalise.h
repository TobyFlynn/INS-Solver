inline void p_multigrid_vec_normalise(const int *order, const double *in,
                                      const double *norm, double *out) {
  const int dg_np = DG_CONSTANTS[(*order - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    out[i] = in[i] / *norm;
  }
}
