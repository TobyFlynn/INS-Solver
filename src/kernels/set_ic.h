inline void set_ic(double *q0, double *q1, double *exQ0,
                   double *exQ1, double *dPdN0, double *dPdN1, double *pRHSex, double *d, double *n0, double *n1) {
  for(int i = 0; i < 15; i++) {
    q0[i] = ic_u;
    q1[i] = ic_v;
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
    dPdN0[i] = 0.0;
    dPdN1[i] = 0.0;
    pRHSex[i] = 0.0;
    d[i] = 0.0;
    n0[i] = 0.0;
    n1[i] = 0.0;
  }
}
