// Box = {bottom left x, bottom left y, upper right x, upper right y}
inline void measure_lift_drag(const DG_FP *box, DG_FP *lift_coeff, DG_FP *drag_coeff,
                              const int *edgeNum, const DG_FP *nx, const DG_FP *ny, 
                              const DG_FP *sJ, const DG_FP *x, const DG_FP *y, 
                              const DG_FP *pr, const DG_FP *ux, const DG_FP *uy, 
                              const DG_FP *vx, const DG_FP *vy) {
  const int edge = *edgeNum;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edge * DG_NPF];

  bool within_box = true;
  for(int i = 0; i < DG_NPF; i++) {
    int fmask_ind = fmask[i];
    if(x[fmask_ind] < box[0] || x[fmask_ind] > box[2]
      || y[fmask_ind] < box[1] || y[fmask_ind] > box[3])
      within_box = false;
  }

  if(!within_box)
    return;

  // Sum each column of emat
  const DG_FP *emat = dg_Emat_kernel + (DG_ORDER - 1) * DG_NUM_FACES * DG_NPF * DG_NP;
  DG_FP w[DG_NPF];
  for(int i = 0; i < DG_NPF; i++) {
    w[i] = 0.0;
    for(int j = 0; j < DG_NPF; j++) {
      int fmask_ind = fmask[j];
      int ind = DG_MAT_IND(fmask_ind, edge * DG_NPF + i, DG_NP, DG_NUM_FACES * DG_NPF);
      w[i] += emat[ind];
    }
  }

  const DG_FP _nx = *nx;
  const DG_FP _ny = *ny;
  const DG_FP _sJ = *sJ;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    DG_FP h_force = -pr[fmask_ind] * _nx + (_nx * 2.0 * ux[fmask_ind] + _ny * (vx[fmask_ind] + uy[fmask_ind])) / r_ynolds;
    DG_FP v_force = -pr[fmask_ind] * _ny + (_nx * (uy[fmask_ind] + vx[fmask_ind]) + _ny * 2.0 * vy[fmask_ind]) / r_ynolds;

    *drag_coeff += -w[i] * _sJ * h_force;
    *lift_coeff += -w[i] * _sJ * v_force;
  }
}
