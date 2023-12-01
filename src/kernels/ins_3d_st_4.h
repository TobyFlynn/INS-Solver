inline void ins_3d_st_4(DG_FP *nx, DG_FP *ny, DG_FP *nz) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP mag = sqrt(nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i]);
    if(mag > 1e-8) {
      nx[i] /= mag;
      ny[i] /= mag;
      nz[i] /= mag;
    }
  }
}
