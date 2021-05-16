inline void glb_ind_kernel(const int **glb, int *glbL, int *glbR) {
  glbL[0] = glb[0][0];
  glbR[0] = glb[1][0];
}
