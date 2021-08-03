inline void poisson_op1(const double *J, const double *Dx, const double *Dy,
                        const double *factor, double *op){
  double tmpX[DG_CUB_NP * DG_NP];
  double tmpY[DG_CUB_NP * DG_NP];

  for(int m = 0; m < DG_CUB_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      tmpX[ind] = J[m] * cubW_g[m] * Dx[ind] * factor[m];
      tmpY[ind] = J[m] * cubW_g[m] * Dy[ind] * factor[m];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      op[c_ind] = 0.0;
      for(int k = 0; k < DG_CUB_NP; k++) {
        // tmpX and tmpY
        int b_ind = k * DG_NP + j;
        // Transpose of Dx and Dy
        int a_ind = k * DG_NP + i;
        op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
      }
    }
  }
}
