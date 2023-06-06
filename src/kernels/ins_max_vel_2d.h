inline void ins_max_vel_2d(DG_FP *max_vel, const DG_FP *u, const DG_FP *v) {
  for(int i = 0; i < DG_NP; i++) {
    DG_FP tmp = u[i] * u[i] + v[i] * v[i];
    if(tmp > *max_vel)
      *max_vel = tmp;
  }
}
