inline void ins_pressure_rhs_2d(const DG_FP *b0, const DG_FP *b1,
                                const DG_FP *dt, const DG_FP *dPdN,
                                DG_FP *dPdNOld, DG_FP *divVelT) {
  DG_FP factor = 1.0 / (*dt);
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor;
  }

  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    dPdNOld[i] = (*b0) * dPdN[i] + (*b1) * dPdNOld[i];
  }
}
