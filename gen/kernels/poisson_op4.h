inline void poisson_op4(const double *mm, const double *factor, double *op, double *tmp){
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      int c_ind = i * 3 + j;
      int mm_ind = j * 3 + i;
      op[c_ind] += mm[mm_ind] * factor[j];
      tmp[mm_ind] = mm[mm_ind];
    }
  }
}
