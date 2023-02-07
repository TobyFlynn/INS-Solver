inline void ins_vis_copy_2d(const double *u, const double *v, double *u_tmp, double *v_tmp) {
  for(int i = 0; i < DG_NP; i++) {
    u_tmp[i] = u[i];
    v_tmp[i] = v[i];
  }
}
