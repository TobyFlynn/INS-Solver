inline void ins_vis_mm_fact_2d(const DG_FP *factor, const DG_FP *rho,
                               DG_FP *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = *factor * rho[i];
  }
}
