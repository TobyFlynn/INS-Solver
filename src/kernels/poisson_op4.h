inline void poisson_op4(const double *mm, const double *factor, double *op, double *tmp){
  for(int i = 0; i < DG_NP; i++) {
    for(int j = 0; j < DG_NP; j++) {
      int c_ind = i * DG_NP + j;
      int mm_ind = j * DG_NP + i;
      op[c_ind] += mm[mm_ind] * factor[j];
      tmp[mm_ind] = mm[mm_ind];
    }
  }
}
