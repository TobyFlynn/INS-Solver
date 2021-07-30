inline void poisson_op4(const double *mm, const double *factor, double *op, double *tmp){
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      int c_ind = i * 6 + j;
      int mm_ind = j * 6 + i;
      op[c_ind] += mm[mm_ind] * factor[j];
      tmp[mm_ind] = mm[mm_ind];
    }
  }
}
