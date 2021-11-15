inline void set_ic(const int *p, double *q0, double *q1) {
  // Get constants for this element's order
  const int dg_np  = DG_CONSTANTS[(*p - 1) * 5];

  for(int i = 0; i < dg_np; i++) {
    q0[i] = ic_u;
    q1[i] = ic_v;
  }
}
