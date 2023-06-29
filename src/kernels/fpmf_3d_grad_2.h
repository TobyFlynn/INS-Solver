inline void fpmf_3d_grad_2(const int *p, const DG_FP *geof,
                           const DG_FP * __restrict__ fact,
                           DG_FP * __restrict__ ux, DG_FP * __restrict__ uy,
                           DG_FP *__restrict__ uz) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  #pragma omp simd
  for(int m = 0; m < dg_np; m++) {
    const DG_FP r = ux[m];
    const DG_FP s = uy[m];
    const DG_FP t = uz[m];
    ux[m] = fact[m] * (geof[RX_IND] * r + geof[SX_IND] * s + geof[TX_IND] * t);
    uy[m] = fact[m] * (geof[RY_IND] * r + geof[SY_IND] * s + geof[TY_IND] * t);
    uz[m] = fact[m] * (geof[RZ_IND] * r + geof[SZ_IND] * s + geof[TZ_IND] * t);
  }
}
