inline void ins_3d_enstrophy_1(DG_FP *enstropy, const DG_FP *in) {
  DG_FP tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    tmp += in[i];
  }
  *enstropy += tmp;
}
