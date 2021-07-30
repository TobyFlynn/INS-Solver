inline void ls_copy(const double *dsdx, const double *dsdy,
                    double *dpldx, double *dprdx, double *dpldy, double *dprdy) {
  for(int i = 0; i < 6; i++) {
    dpldx[i] = dsdx[i];
    dprdx[i] = dsdx[i];
    dpldy[i] = dsdy[i];
    dprdy[i] = dsdy[i];
  }
}
