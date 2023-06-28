inline void ins_3d_enstrophy_0(const DG_FP *curl0, const DG_FP *curl1,
                               DG_FP *curl2) {
  for(int i = 0; i < DG_NP; i++) {
    curl2[i] = curl0[i] * curl0[i] + curl1[i] * curl1[i] + curl2[i] * curl2[i];
  }
}
