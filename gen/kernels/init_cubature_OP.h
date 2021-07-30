inline void init_cubature_OP(const double *J, const double *Dx,
                             const double *Dy, double *temp, double *temp2) {
  for(int m = 0; m < 16; m++) {
    for(int n = 0; n < 6; n++) {
      int ind = m * 6 + n;
      temp[ind] = J[m] * cubW_g[m] * Dx[ind];
      temp2[ind] = J[m] * cubW_g[m] * Dy[ind];
    }
  }
}
