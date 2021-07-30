inline void poisson_op1(const double *J, const double *Dx, const double *Dy,
                        const double *factor, double *op){
  double tmpX[36 * 10];
  double tmpY[36 * 10];

  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      tmpX[ind] = J[m] * cubW_g[m] * Dx[ind] * factor[m];
      tmpY[ind] = J[m] * cubW_g[m] * Dy[ind] * factor[m];
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op[c_ind] = 0.0;
      for(int k = 0; k < 36; k++) {
        // tmpX and tmpY
        int b_ind = k * 10 + j;
        // Transpose of Dx and Dy
        int a_ind = k * 10 + i;
        op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
      }
    }
  }
}
