inline void cg3_update_p(const DG_FP *beta, const DG_FP *r0, const DG_FP *r1,
                         const DG_FP *r2, DG_FP *p0, DG_FP *p1, DG_FP *p2) {
  for(int i = 0; i < DG_NP; i++) {
    p0[i] = r0[i] + *beta * p0[i];
    p1[i] = r1[i] + *beta * p1[i];
    p2[i] = r2[i] + *beta * p2[i];
  }
}