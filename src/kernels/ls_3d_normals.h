inline void ls_3d_normals(const DG_FP *s, DG_FP *nx, DG_FP *ny, DG_FP *nz) {
  for(int i = 0; i < DG_NP; i++) {
    DG_FP size = sqrt(nx[i] * nx[i] + ny[i] * ny[i] + nz[i] * nz[i]);
    nx[i] /= size;
    ny[i] /= size;
    nz[i] /= size;
    // if(s[i] < 0.0) {
    //   nx[i] = -nx[i];
    //   ny[i] = -ny[i];
    // }
  }
}
