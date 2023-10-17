inline void diff_ex(const DG_FP *x, const DG_FP *y, DG_FP *val, DG_FP *vis) {
  for(int i = 0; i < DG_NP; i++) {
    if(y[i] > 0.0)
      vis[i] = fmax(0.1, 1.0 - 0.1 * y[i]);
    else
      vis[i] = fmax(0.1, 1.0 + y[i]);
    val[i] = x[i] * x[i] + y[i] * y[i] - 10.0;
    val[i] = fmax(0.0, -val[i]);
  }
}