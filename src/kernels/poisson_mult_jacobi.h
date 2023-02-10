inline void poisson_mult_jacobi(const int *p, const double *op, double *rhs) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    const int op_ind = i * dg_np + i;
    rhs[i] = rhs[i] / op[op_ind];
  }
}
