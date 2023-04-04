inline void ins_3d_proj_pen(const DG_FP *factor, const DG_FP *u, const DG_FP *v,
                            const DG_FP *w, const DG_FP *h, DG_FP *pen) {
  DG_FP max_vel = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP vel = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
    if(vel > max_vel)
      max_vel = vel;
  }
  max_vel = sqrt(max_vel);

  *pen = *factor * (*h) * max_vel;
}
