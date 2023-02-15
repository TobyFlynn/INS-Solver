inline void poisson_2d_apply_bc(const int *p, const int *bedgeNum,
                                const DG_FP *op, const DG_FP *bc, DG_FP *rhs) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  int exInd = *bedgeNum * DG_GF_NP;

  for(int m = 0; m < dg_np; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      // int ind = m + n * dg_np;
      int ind = DG_MAT_IND(m, n, dg_np, DG_GF_NP);
      rhs[m] += op[ind] * bc[exInd + n];
    }
  }
}
