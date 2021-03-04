inline void init_cubature_OP(const double *J, const double *Dx, const double *Dy,
                             double *temp, double *temp2) {
  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      temp[ind] = J[m] * cubW[m] * Dx[ind];
      temp2[ind] = J[m] * cubW[m] * Dy[ind];
    }
  }
}
