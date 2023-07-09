inline void fpmf_2d_grad(const int *p, const DG_FP *geof,
                         const DG_FP * __restrict__ fact,
                         DG_FP * __restrict__ ux, DG_FP * __restrict__ uy) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  for(int m = 0; m < dg_np; m++) {
    const DG_FP r = ux[m];
    const DG_FP s = uy[m];
    ux[m] = fact[m] * (geof[RX_IND] * r + geof[SX_IND] * s);
    uy[m] = fact[m] * (geof[RY_IND] * r + geof[SY_IND] * s);
  }
}
