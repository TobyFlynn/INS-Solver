inline void constant_mult(const DG_FP *constant, DG_FP *val) {
  for(int i = 0; i < DG_NP; i++) {
    val[i] *= *constant;
  }
}