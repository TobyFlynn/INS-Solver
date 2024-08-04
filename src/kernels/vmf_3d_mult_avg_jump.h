inline void vmf_3d_mult_avg_jump(const int *u_bc_type, const int *v_bc_type, const int *w_bc_type, 
                    const int *faceNum, const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, 
                    const DG_FP *u, const DG_FP *u_x, const DG_FP *u_y, const DG_FP *u_z, 
                    const DG_FP *v, const DG_FP *v_x, const DG_FP *v_y, const DG_FP *v_z,
                    const DG_FP *w, const DG_FP *w_x, const DG_FP *w_y, const DG_FP *w_z,
                    DG_FP *u_jump, DG_FP *u_x_avg, DG_FP *u_y_avg, DG_FP *u_z_avg, DG_FP *v_jump, 
                    DG_FP *v_x_avg, DG_FP *v_y_avg, DG_FP *v_z_avg, DG_FP *w_jump, DG_FP *w_x_avg, 
                    DG_FP *w_y_avg, DG_FP *w_z_avg) {
  if(*u_bc_type == 2 && *v_bc_type == 2 && *w_bc_type == 2) {
    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskInd = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * DG_NPF + j];
      const int ind = faceNum[0] * DG_NPF + j;

      const DG_FP projected_u = *nx * *nx * u[fmaskInd] + *nx * *ny * v[fmaskInd] + *nx * *nz * w[fmaskInd];
      const DG_FP projected_v = *nx * *ny * u[fmaskInd] + *ny * *ny * v[fmaskInd] + *ny * *nz * w[fmaskInd];
      const DG_FP projected_w = *nx * *nz * u[fmaskInd] + *ny * *nz * v[fmaskInd] + *nz * *nz * w[fmaskInd];
      const DG_FP projected_ux = *nx * *nx * u_x[fmaskInd] + *nx * *ny * v_x[fmaskInd] + *nx * *nz * w_x[fmaskInd];
      const DG_FP projected_vx = *nx * *ny * u_x[fmaskInd] + *ny * *ny * v_x[fmaskInd] + *ny * *nz * w_x[fmaskInd];
      const DG_FP projected_wx = *nx * *nz * u_x[fmaskInd] + *ny * *nz * v_x[fmaskInd] + *nz * *nz * w_x[fmaskInd];
      const DG_FP projected_uy = *nx * *nx * u_y[fmaskInd] + *nx * *ny * v_y[fmaskInd] + *nx * *nz * w_y[fmaskInd];
      const DG_FP projected_vy = *nx * *ny * u_y[fmaskInd] + *ny * *ny * v_y[fmaskInd] + *ny * *nz * w_y[fmaskInd];
      const DG_FP projected_wy = *nx * *nz * u_y[fmaskInd] + *ny * *nz * v_y[fmaskInd] + *nz * *nz * w_y[fmaskInd];
      const DG_FP projected_uz = *nx * *nx * u_z[fmaskInd] + *nx * *ny * v_z[fmaskInd] + *nx * *nz * w_z[fmaskInd];
      const DG_FP projected_vz = *nx * *ny * u_z[fmaskInd] + *ny * *ny * v_z[fmaskInd] + *ny * *nz * w_z[fmaskInd];
      const DG_FP projected_wz = *nx * *nz * u_z[fmaskInd] + *ny * *nz * v_z[fmaskInd] + *nz * *nz * w_z[fmaskInd];
      
      u_jump[ind] += 2.0 * projected_u;
      v_jump[ind] += 2.0 * projected_v;
      w_jump[ind] += 2.0 * projected_w;
      u_x_avg[ind] += projected_ux;
      v_x_avg[ind] += projected_vx;
      w_x_avg[ind] += projected_wx;
      u_y_avg[ind] += projected_uy;
      v_y_avg[ind] += projected_vy;
      w_y_avg[ind] += projected_wy;
      u_z_avg[ind] += projected_uz;
      v_z_avg[ind] += projected_vz;
      w_z_avg[ind] += projected_wz;
    }
    return;
  }
  
  if(*u_bc_type == 1) {
    // Do nothing
  } else if(*u_bc_type == 0) {
    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskInd = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * DG_NPF + j];
      const int ind = faceNum[0] * DG_NPF + j;
      u_jump[ind] += 2.0 * u[fmaskInd];
      u_x_avg[ind] += u_x[fmaskInd];
      u_y_avg[ind] += u_y[fmaskInd];
      u_z_avg[ind] += u_z[fmaskInd];
    }
  }

  if(*v_bc_type == 1) {
    // Do nothing
  } else if(*v_bc_type == 0) {
    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskInd = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * DG_NPF + j];
      const int ind = faceNum[0] * DG_NPF + j;
      v_jump[ind] += 2.0 * v[fmaskInd];
      v_x_avg[ind] += v_x[fmaskInd];
      v_y_avg[ind] += v_y[fmaskInd];
      v_z_avg[ind] += v_z[fmaskInd];
    }
  }

  if(*w_bc_type == 1) {
    // Do nothing
  } else if(*w_bc_type == 0) {
    for(int j = 0; j < DG_NPF; j++) {
      const int fmaskInd = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + faceNum[0] * DG_NPF + j];
      const int ind = faceNum[0] * DG_NPF + j;
      w_jump[ind] += 2.0 * w[fmaskInd];
      w_x_avg[ind] += w_x[fmaskInd];
      w_y_avg[ind] += w_y[fmaskInd];
      w_z_avg[ind] += w_z[fmaskInd];
    }
  }
}