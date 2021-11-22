inline void ls_add_diff(const double *diff, double *rk) {
  for(int i = 0; i < DG_NP; i++) {
    rk[i] = rk[i] + diff[i];
  }
}
