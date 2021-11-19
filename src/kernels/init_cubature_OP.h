inline void init_cubature_OP(const double *J, const double *Dx, const double *Dy,
                             double *temp, double *temp2) {
  for(int m = 0; m < DG_CUB_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      temp[ind] = J[m] * cubW_g[m] * Dx[ind];
      temp2[ind] = J[m] * cubW_g[m] * Dy[ind];
    }
  }
}
