inline void ls_add_diff(const double *diff, double *rk, double *dsldx,
                        double *dsrdx, double *dsldy, double *dsrdy) {
  for(int i = 0; i < 3; i++) {
    rk[i] = rk[i] + diff[i];
  }

  for(int i = 0; i < 9; i++) {
    dsldx[i] = 0.0;
    dsrdx[i] = 0.0;
    dsldy[i] = 0.0;
    dsrdy[i] = 0.0;
  }
}
