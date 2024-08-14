inline void euler_2d_vol(const DG_FP *q0, const DG_FP *q1, const DG_FP *q2, 
                         const DG_FP *q3, DG_FP *f0, DG_FP *f1, DG_FP *f2, 
                         DG_FP *f3, DG_FP *g0, DG_FP *g1, DG_FP *g2, 
                         DG_FP *g3) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP rho = q0[i];
    const DG_FP rhou = q1[i];
    const DG_FP rhov = q2[i];
    const DG_FP energy = q3[i];
    const DG_FP u = rhou / rho;
    const DG_FP v = rhov / rho;
    const DG_FP p = (1.4 - 1.0) * (energy - 0.5 * (rhou * u + rhov * v));

    f0[i] = rhou;
    f1[i] = rhou * u + p;
    f2[i] = rhou * v;
    f3[i] = u * (energy + p);

    g0[i] = rhov;
    g1[i] = rhou * v;
    g2[i] = rhov * v + p;
    g3[i] = v * (energy + p);
  }
}