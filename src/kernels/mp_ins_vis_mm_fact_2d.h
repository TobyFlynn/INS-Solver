inline void mp_ins_vis_mm_fact_2d(const double *factor, const double *rho,
                                  double *out) {
  for(int i = 0; i < DG_NP; i++) {
    out[i] = *factor * rho[i];
  }
}
