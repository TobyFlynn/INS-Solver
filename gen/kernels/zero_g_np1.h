inline void zero_g_np1(double *g_np0) {
  for(int i = 0; i < 18; i++) {
    g_np0[i] = 0.0;
  }
}
