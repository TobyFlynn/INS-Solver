inline void zero_cub_surf_2d(DG_FP *val) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    val[i] = 0.0;
  }
}
