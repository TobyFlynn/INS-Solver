inline void ins_3d_vis_1(const DG_FP *t, const int *bc_type, const int *faceNum,
                         const DG_FP *x, const DG_FP *y, const DG_FP *z,
                         int *type, DG_FP *bcs) {
  if(*bc_type == 1)
    *type = 1;
  else
    *type = 0;

  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const DG_FP PI = 3.141592653589793238463;
  if(*bc_type == 0) {
    for(int i = 0; i < DG_NPF; i++) {
      // bcs[i] = sin(PI * (*t)) * (y[fmask[i]] * (1.0 - y[fmask[i]]))  * (z[fmask[i]] * (1.0 - z[fmask[i]]));
      bcs[i] = sin(PI * (*t));
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      bcs[i] = 0.0;
    }
  }
}
