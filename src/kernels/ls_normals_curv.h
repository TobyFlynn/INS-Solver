inline void ls_normals_curv(const DG_FP *alpha, const DG_FP *s, const DG_FP *x, const DG_FP *y, DG_FP *nx, DG_FP *ny, DG_FP *curv) {
  for(int i = 0; i < DG_NP; i++) {
    // nx[i] = x[i];
    // ny[i] = y[i];
    // DG_FP size = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    // if(size > 1e-5) {
    //   nx[i] /= size;
    //   ny[i] /= size;
    // }
    curv[i] = 1.0 / 7.5;
  }
}
