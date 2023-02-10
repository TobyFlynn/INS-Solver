inline void ls_update_order(const int *p, const double *alpha, const double *s,
                            int *order) {
  // Get constants for this element's order
  const int dg_np = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];

  *order = DG_ORDER;
  for(int i = 0; i < dg_np; i++) {
    if(fabs(s[i]) < *alpha) {
      // *order = 1;
    }
  }
}
