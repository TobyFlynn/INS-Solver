inline void poisson_pre(const double *in, const double *pre, double *out) {
  /*
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      int mm_ind = j * 15 + i;
      op[c_ind] += mm[mm_ind] * factor[j];
    }
  }

  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    rhs[m] = 0.0;
    for(int n = 0; n < 15; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
  */

  for(int i = 0; i < 15; i++) {
    out[i] = 0.0;
    // out[i] = in[i];
    for(int j = 0; j < 15; j++) {
      int mm_ind = j * 15 + i;
      // out[i] += mmInv[mm_ind] * (1.0 / factor[j]) * in[j];
      out[i] += pre[mm_ind] * in[j];
    }
  }
}
