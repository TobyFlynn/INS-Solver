inline void poisson_op1(const double *J, const double *Dx, const double *Dy,
                        const double *factor, double *op){
  double tmpX[46*15];
  double tmpY[46*15];

  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      tmpX[ind] = J[m] * cubW_g[m] * Dx[ind] * factor[m];
      tmpY[ind] = J[m] * cubW_g[m] * Dy[ind] * factor[m];
    }
  }

  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op[c_ind] = 0.0;
      for(int k = 0; k < 46; k++) {
        // tmpX and tmpY
        int b_ind = k * 15 + j;
        // Transpose of Dx and Dy
        int a_ind = k * 15 + i;
        op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
      }
    }
  }
}
