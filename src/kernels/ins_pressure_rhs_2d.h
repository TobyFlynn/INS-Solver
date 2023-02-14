inline void ins_pressure_rhs_2d(const DG_FP *b0, const DG_FP *b1,
                                const DG_FP *dt, const int *p, const DG_FP *J,
                                const DG_FP *sJ, const DG_FP *dPdN,
                                DG_FP *dPdNOld, DG_FP *divVelT) {
  DG_FP factor = 1.0 / (*dt);
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor * J[i];
  }

  const DG_FP *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];
  for(int i = 0; i < DG_G_NP; i++) {
    dPdNOld[i] = gaussW[i % DG_GF_NP] * sJ[i] * ((*b0) * dPdN[i] + (*b1) * dPdNOld[i]);
  }
}
