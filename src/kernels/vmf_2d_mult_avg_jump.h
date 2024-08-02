inline void vmf_2d_mult_avg_jump(const int *order, const int *u_bc_type,
                                 const int *v_bc_type, const int *faceNum,
                                 const DG_FP *nx, const DG_FP *ny, 
                                 const DG_FP *u, const DG_FP *u_x, const DG_FP *u_y,
                                 const DG_FP *v, const DG_FP *v_x, const DG_FP *v_y, 
                                 DG_FP *u_jump, DG_FP *u_x_avg, DG_FP *u_y_avg,
                                 DG_FP *v_jump, DG_FP *v_x_avg, DG_FP *v_y_avg) {
  if(*u_bc_type == 2 && *v_bc_type == 2) {
    const int p = order[0];
    const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

    for(int j = 0; j < dg_npf; j++) {
      const int fmaskInd = FMASK[(p - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * dg_npf + j];
      const int ind = faceNum[0] * dg_npf + j;

      const DG_FP projected_u = *nx * *nx * u[fmaskInd] + *nx * *ny * v[fmaskInd];
      const DG_FP projected_v = *nx * *ny * u[fmaskInd] + *ny * *ny * v[fmaskInd];
      const DG_FP projected_ux = *nx * *nx * u_x[fmaskInd] + *nx * *ny * v_x[fmaskInd];
      const DG_FP projected_vx = *nx * *ny * u_x[fmaskInd] + *ny * *ny * v_x[fmaskInd];
      const DG_FP projected_uy = *nx * *nx * u_y[fmaskInd] + *nx * *ny * v_y[fmaskInd];
      const DG_FP projected_vy = *nx * *ny * u_y[fmaskInd] + *ny * *ny * v_y[fmaskInd];
      
      u_jump[ind] += 2.0 * projected_u;
      v_jump[ind] += 2.0 * projected_v;
      u_x_avg[ind] += projected_ux;
      v_x_avg[ind] += projected_vx;
      u_y_avg[ind] += projected_uy;
      v_y_avg[ind] += projected_vy;
    }
    return;
  }
  
  if(*u_bc_type == 1) {
    // Do nothing
  } else if(*u_bc_type == 0) {
    const int p = order[0];
    const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

    for(int j = 0; j < dg_npf; j++) {
      const int fmaskInd = FMASK[(p - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * dg_npf + j];
      const int ind = faceNum[0] * dg_npf + j;
      u_jump[ind] += 2.0 * u[fmaskInd];
      u_x_avg[ind] += u_x[fmaskInd];
      u_y_avg[ind] += u_y[fmaskInd];
    }
  }

  if(*v_bc_type == 1) {
    // Do nothing
  } else if(*v_bc_type == 0) {
    const int p = order[0];
    const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

    for(int j = 0; j < dg_npf; j++) {
      const int fmaskInd = FMASK[(p - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * dg_npf + j];
      const int ind = faceNum[0] * dg_npf + j;
      v_jump[ind] += 2.0 * v[fmaskInd];
      v_x_avg[ind] += v_x[fmaskInd];
      v_y_avg[ind] += v_y[fmaskInd];
    }
  }
}