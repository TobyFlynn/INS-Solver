inline void measure_min_max_interface(DG_FP *min_x, DG_FP *min_y, DG_FP *max_x, 
                                      DG_FP *max_y, const DG_FP *x, const DG_FP *y, 
                                      const DG_FP *ls) {
  for(int i = 0; i < DG_NP; i++) {
    if(ls[i] < 0.0) {
      *min_x = fmin(*min_x, x[i]);
      *min_y = fmin(*min_y, y[i]);
      *max_x = fmax(*max_x, x[i]);
      *max_y = fmax(*max_y, y[i]);
    }
  }
}