inline void mp_ins_2d_vis_0(const DG_FP *fact0, const DG_FP *fact1,
                            const DG_FP *rho, DG_FP *rhs0, DG_FP *rhs1,
                            DG_FP *mm_factor) {
  for(int i = 0; i < DG_NP; i++) {
    rhs0[i] *= *fact0 * rho[i];
    rhs1[i] *= *fact0 * rho[i];
    mm_factor[i] = *fact1 * rho[i];
  }
}
