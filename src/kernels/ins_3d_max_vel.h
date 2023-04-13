inline void ins_3d_max_vel(DG_FP *max_vel, const DG_FP *u, const DG_FP *v, const DG_FP *w) {
  for(int i = 0; i < DG_NP; i++) {
    DG_FP tmp = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
    if(tmp > *max_vel)
      *max_vel = tmp;
  }
}
