inline void euler_2d_vol(const DG_FP *geof, DG_FP *f0, DG_FP *f1, 
                         DG_FP *f2, DG_FP *f3, DG_FP *g0, DG_FP *g1, 
                         DG_FP *g2, DG_FP *g3) {
  const DG_FP rx = geof[RX_IND];
  const DG_FP sx = geof[SX_IND];
  const DG_FP ry = geof[RY_IND];
  const DG_FP sy = geof[SY_IND];
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    const DG_FP rho = f0[i];
    const DG_FP rhou = f1[i];
    const DG_FP rhov = f2[i];
    const DG_FP energy = f3[i];
    const DG_FP u = rhou / rho;
    const DG_FP v = rhov / rho;
    const DG_FP p = (1.4 - 1.0) * (energy - 0.5 * (rhou * u + rhov * v));

    const DG_FP _f0 = rhou;
    const DG_FP _f1 = rhou * u + p;
    const DG_FP _f2 = rhou * v;
    const DG_FP _f3 = u * (energy + p);

    const DG_FP _g0 = rhov;
    const DG_FP _g1 = rhou * v;
    const DG_FP _g2 = rhov * v + p;
    const DG_FP _g3 = v * (energy + p);

    // F used for dr, G used for ds
    f0[i] = rx * _f0 + ry * _g0;
    g0[i] = sx * _f0 + sy * _g0;
    f1[i] = rx * _f1 + ry * _g1;
    g1[i] = sx * _f1 + sy * _g1;
    f2[i] = rx * _f2 + ry * _g2;
    g2[i] = sx * _f2 + sy * _g2;
    f3[i] = rx * _f3 + ry * _g3;
    g3[i] = sx * _f3 + sy * _g3;
  }
}