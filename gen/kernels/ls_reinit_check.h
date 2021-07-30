inline void ls_reinit_check(const double *alpha, const double *s,
                            const double *dsdx, const double *dsdy,
                            double *res, int *count) {
  for(int i = 0; i < 10; i++) {
    if(fabs(s[i]) < (*alpha)) {
      *res += dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
      *count += 1;
    }
  }
}
