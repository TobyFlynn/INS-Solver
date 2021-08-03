inline void init_cubature_OP(const double *J, const double *Dx,
                             const double *Dy, double *temp, double *temp2) {
  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      temp[ind] = J[m] * cubW_g[m] * Dx[ind];
      temp2[ind] = J[m] * cubW_g[m] * Dy[ind];
    }
  }
}
