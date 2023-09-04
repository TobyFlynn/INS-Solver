inline void ins_2d_st_5(DG_FP *nx, DG_FP *ny) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP mag = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    if(mag > 1e-8) {
      nx[i] /= mag;
      ny[i] /= mag;
    } 
  }
}