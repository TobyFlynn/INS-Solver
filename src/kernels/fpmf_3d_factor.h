inline void fpmf_3d_factor(const int *p, const DG_FP * __restrict__ fact,
                           DG_FP * __restrict__ ux, DG_FP * __restrict__ uy,
                           DG_FP *__restrict__ uz) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  #pragma omp simd
  for(int m = 0; m < dg_np; m++) {
    ux[m] *= fact[m];
    uy[m] *= fact[m];
    uz[m] *= fact[m];
  }
}
