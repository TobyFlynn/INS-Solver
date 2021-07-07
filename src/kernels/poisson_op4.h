inline void poisson_op4(const double *mm, const double *factor, double *op){
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      int mm_ind = j * 15 + i;
      op[c_ind] += mm[mm_ind] * factor[j];
    }
  }
}
