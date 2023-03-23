inline void petsc_pre_jacobi(const DG_FP *diag, const DG_FP *in, DG_FP* out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = in[i] / diag[i];
  }
}