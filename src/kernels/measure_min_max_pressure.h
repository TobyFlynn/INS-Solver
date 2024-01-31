inline void measure_min_max_pressure(const DG_FP *alpha, DG_FP *min_pr, DG_FP *max_pr, const DG_FP *pr, const DG_FP *ls) {
  for(int i = 0; i < DG_NP; i++) {
    if(ls[i] < fabs(*alpha)) {
      *min_pr = fmin(*min_pr, pr[i]);
      *max_pr = fmax(*max_pr, pr[i]);
    }
  }
}