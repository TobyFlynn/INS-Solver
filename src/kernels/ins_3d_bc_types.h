inline void ins_3d_bc_types(const int *faceNum, const DG_FP *x, const DG_FP *y, 
                            const DG_FP *z, int *type) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  bool x_same = true;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind_0 = fmask[0];
    const int fmask_ind_i = fmask[i];
    if(!fp_equal(x[fmask_ind_0], x[fmask_ind_i]))
      x_same = false;
  }
  bool x_less_than_zero = true;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind_i = fmask[i];
    if(x[fmask_ind_i] > 1e-8)
      x_less_than_zero = false;
  }
  bool x_near_inflow = true;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind_i = fmask[i];
    if(x[fmask_ind_i] > -LW_INLET_LENGTH + 0.1)
      x_near_inflow = false;
  }
  bool x_near_outflow = true;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind_i = fmask[i];
    if(x[fmask_ind_i] < LW_LENGTH_SHORT - 0.1)
      x_near_outflow = false;
  }
  bool face_on_outer_bc = false;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind_i = fmask[i];
    if(y[fmask_ind_i] * y[fmask_ind_i] + z[fmask_ind_i] * z[fmask_ind_i] > LW_RADIUS * LW_RADIUS - 5e-1)
      face_on_outer_bc = true;
  }
  
  if(x_same && x_near_inflow) {
    *type = LW_INFLOW_BC;
  } else if(x_less_than_zero) {
    *type = LW_NO_SLIP_WALL_BC;
  } else if(x_same && x_near_outflow) {
    *type = LW_OUTFLOW_BC;
  } else if(!face_on_outer_bc) {
    *type = LW_NO_SLIP_WALL_BC;
  } else {
    *type = LW_SLIP_WALL_BC;
  }
}
