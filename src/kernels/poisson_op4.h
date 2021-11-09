inline void poisson_op4(const double *cJ, const double *factor, double *op, double *tmp){
  double cTmp[DG_CUB_NP * DG_NP];
  double mm[DG_NP * DG_NP];
  for(int m = 0; m < DG_CUB_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_CUB_NP;
      cTmp[ind] = factor[m] * cJ[m] * cubW_g[m] * cubV_g[ind];
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i + j * DG_NP;
      mm[c_ind] = 0.0;
      for(int k = 0; k < DG_CUB_NP; k++) {
        int b_ind = k + j * DG_CUB_NP;
        // Transpose
        int a_ind = k + i * DG_CUB_NP;

        mm[c_ind] += cubV_g[a_ind] * cTmp[b_ind];
      }
    }
  }

  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i + j * DG_NP;
      op[c_ind] += mm[c_ind];
      tmp[c_ind] = mm[c_ind];
      // int mm_ind = j * DG_NP + i;
      // op[c_ind] += mm[c_ind] * factor[j];
      // tmp[c_ind] = mm[c_ind] * factor[j];
    }
  }
}
