inline void poisson_op4(const double *cJ, const double *factor, double *op, double *tmp){
  double cTmp[36 * 10];
  double mm[10 * 10];
  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      cTmp[ind] = factor[m] * cJ[m] * cubW_g[m] * cubV_g[ind];
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      mm[c_ind] = 0.0;
      for(int k = 0; k < 36; k++) {
        int b_ind = k * 10 + j;
        // Transpose
        int ind = i * 36 + k;
        int a_ind = ((ind * 10) % (10 * 36)) + (ind / 36);

        mm[c_ind] += cubV_g[b_ind] * cTmp[a_ind];
      }
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op[c_ind] += mm[c_ind];
      tmp[c_ind] = mm[c_ind];
      // int mm_ind = j * 10 + i;
      // op[c_ind] += mm[c_ind] * factor[j];
      // tmp[c_ind] = mm[c_ind] * factor[j];
    }
  }
}
