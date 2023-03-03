inline void euler_2d_flux(const DG_FP *q0, const DG_FP *q1, 
                          const DG_FP *q2, const DG_FP *q3,
                          DG_FP *f0, DG_FP *f1, DG_FP *f2, DG_FP *f3,
                          DG_FP *g0, DG_FP *g1, DG_FP *g2, DG_FP *g3) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP u = q1[i] / q0[i];
    const DG_FP v = q2[i] / q0[i];
    const DG_FP p = (gamma_e - 1.0) * (q3[i] - 0.5 * (q1[i] * u + q2[i] * v));

    f0[i] = q1[i];
    f1[i] = q1[i] * u + p;
    f2[i] = q2[i] * u;
    f3[i] = u * (q3[i] + p);

    g0[i] = q2[i];
    g1[i] = q1[i] * v;
    g2[i] = q2[i] * v + p;
    g3[i] = v * (q3[i] + p);
  }
}