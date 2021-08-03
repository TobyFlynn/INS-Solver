inline void zero_g_np(double *g_np0, double *g_np1) {
  for(int i = 0; i < 18; i++) {
    g_np0[i] = 0.0;
    g_np1[i] = 0.0;
  }
}
