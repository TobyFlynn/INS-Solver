inline void ls_reinit_check(const int *p, const DG_FP *alpha, const DG_FP *s,
                            const DG_FP *dsdx, const DG_FP *dsdy,
                            DG_FP *res, int *count) {
  // Get constants for this element's order
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  for(int i = 0; i < dg_np; i++) {
    if(fabs(s[i]) < (*alpha)) {
      *res += dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
      *count += 1;
    }
  }
}
