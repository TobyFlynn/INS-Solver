inline void poisson_jacobi(const int *p, const double *op, double *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    const int op_ind = i * dg_np + i;
    rhs[i] = rhs[i] / op[op_ind];
  }
}