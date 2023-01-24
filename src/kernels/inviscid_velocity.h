inline void inviscid_velocity(const double *factor, const int *p, const double *uTT, 
                              const double *vTT, double *u, double *v) {
  const int dg_np = DG_CONSTANTS[(*p - 1) * 5];
  for(int i = 0; i < dg_np; i++) {
    u[i] = uTT[i] / *factor;
    v[i] = vTT[i] / *factor;
  }
}