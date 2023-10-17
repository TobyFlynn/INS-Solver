inline void ins_3d_set_ic(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                          DG_FP *u0, DG_FP *v0, DG_FP *w0, DG_FP *u1,
                          DG_FP *v1, DG_FP *w1) {
  for(int i = 0; i < DG_NP; i++) {
    /*
    DG_FP tmp = y[i] * y[i] + z[i] * z[i] - LW_INLET_RADIUS * LW_INLET_RADIUS;
    tmp = fabs(tmp);
    tmp = fmin(1.0, tmp / 0.1);
    u0[i] = tmp * fmax(0.0, 0.0 - (x[i] + 1.0));
    v0[i] = 0.0;
    w0[i] = 0.0;
    */

    u0[i] = sin(x[i]) * cos(y[i]) * cos(z[i]);
    v0[i] = -cos(x[i]) * sin(y[i]) * cos(z[i]);
    w0[i] = 0.0;

    u1[i] = u0[i];
    v1[i] = v0[i];
    w1[i] = w0[i];
  }
}
