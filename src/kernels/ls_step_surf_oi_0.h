inline void ls_step_surf_oi_0(const DG_FP *s, DG_FP *out) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    int fmask_ind = fmask[i];
    out[i] = s[fmask_ind];
  }
}
