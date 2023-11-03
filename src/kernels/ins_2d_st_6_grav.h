inline void ins_2d_st_6_grav(const DG_FP *curv, DG_FP *st_x, DG_FP *st_y) {
  for(int i = 0; i < DG_NP; i++) {
    st_x[i] *= curv[i] / weber;
    st_y[i] *= curv[i] / weber;
    // Apply gravity
    st_y[i] += 1.0 / (froude * froude);
  }
}
