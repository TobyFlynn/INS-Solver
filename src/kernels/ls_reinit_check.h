inline void ls_reinit_check(const DG_FP *alpha, const DG_FP *s,
                            const DG_FP *dsdx, const DG_FP *dsdy,
                            DG_FP *res, int *count) {
  for(int i = 0; i < DG_NP; i++) {
    if(fabs(s[i]) < (*alpha)) {
      *res += dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
      *count += 1;
    }
  }
}
