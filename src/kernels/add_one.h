inline void add_one(DG_FP *dat) {
  for(int i = 0; i < DG_NP; i++) {
    dat[i] += 1.0;
  }
}