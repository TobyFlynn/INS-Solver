inline void fpmf_3d_grad_2(const int *p, const DG_FP *rx, const DG_FP *sx,
                           const DG_FP *tx, const DG_FP *ry, const DG_FP *sy,
                           const DG_FP *ty, const DG_FP *rz, const DG_FP *sz,
                           const DG_FP *tz, const DG_FP * __restrict__ fact,
                           DG_FP * __restrict__ ux, DG_FP * __restrict__ uy,
                           DG_FP *__restrict__ uz) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  #pragma omp simd
  for(int m = 0; m < dg_np; m++) {
    const DG_FP r = ux[m];
    const DG_FP s = uy[m];
    const DG_FP t = uz[m];
    ux[m] = fact[m] * (rx[0] * r + sx[0] * s + tx[0] * t);
    uy[m] = fact[m] * (ry[0] * r + sy[0] * s + ty[0] * t);
    uz[m] = fact[m] * (rz[0] * r + sz[0] * s + tz[0] * t);
  }
}
