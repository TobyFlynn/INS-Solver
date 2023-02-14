inline void ins_vis_copy_2d(const DG_FP *u, const DG_FP *v, DG_FP *u_tmp, DG_FP *v_tmp) {
  for(int i = 0; i < DG_NP; i++) {
    u_tmp[i] = u[i];
    v_tmp[i] = v[i];
  }
}
