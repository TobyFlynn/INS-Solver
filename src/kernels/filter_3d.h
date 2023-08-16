inline void filter_3d(const double *_alpha, const int *_Nc, const int *_sp, DG_FP *modal) {
  const int Nc = *_Nc;
  const int sp = *_sp;
  const DG_FP alpha = (DG_FP)*_alpha;
  int modal_ind = 0;
  for(int i = 0; i < DG_ORDER + 1; i++) {
    for(int j = 0; j < DG_ORDER - i + 1; j++) {
      for(int k = 0; k < DG_ORDER - i - j + 1; k++) {
        int o = i + j + k;
        if(o > Nc)
          modal[modal_ind] *= exp(-alpha*pow((DG_FP)(o - Nc)/(DG_FP)(DG_ORDER-Nc+1),sp));
        modal_ind++;
      }
    }
  }
}