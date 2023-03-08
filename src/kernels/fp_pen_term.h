inline void fp_pen_term(const int *p, const DG_FP *alpha, const DG_FP *h,
                        const DG_FP *factor, const DG_FP *s, DG_FP *pen) {
  const DG_FP PI = 3.141592653589793238463;
  // Get constants
  // Using same Gauss points so should be able to replace dg_gf_npL and
  // dg_gf_npR with DG_GF_NP
  const int dg_np  = DG_CONSTANTS[(p[0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(p[0] - 1) * DG_NUM_CONSTANTS + 1];

  DG_FP max_hinv = 1.0;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP deltaL = 0.0;
    DG_FP minDelta = 0.5 * (PI / *alpha) * (1.0 / (cosh(PI * 1.5) * cosh(PI * 1.5)));
    if(fabs(s[i]) < 1.5 * *alpha) {
      // deltaL = 0.5 * ((1.0 / *alpha) + (1.0 / *alpha) * cos(PI * s[0][indL] / *alpha));
      deltaL = 0.5 * (PI / *alpha) * (1.0 / (cosh(PI * s[i] / *alpha) * cosh(PI * s[i] / *alpha))) - minDelta;
    }

    DG_FP delta = 1.0 + 1e2 * deltaL;
    // DG_FP delta = 1.0;

    pen[i] = 2.0 * max_hinv * delta * (p[0] + 1) * (p[0] + 2) / factor[i];
  }
}
