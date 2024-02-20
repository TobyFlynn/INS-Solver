inline void measure_enstrophy(const DG_FP *w, const DG_FP *geof, const DG_FP *curl0, 
                              const DG_FP *curl1, DG_FP *curl2) {
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    curl2[i] = curl0[i] * curl0[i] + curl1[i] * curl1[i] + curl2[i] * curl2[i];
    curl2[i] *= w[i] * geof[J_IND];
  }
}