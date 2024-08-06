inline void vmf_3d_mult_flux(const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, 
                  const DG_FP *sJ, const DG_FP *tau, DG_FP *u_jump, DG_FP *u_x_avg, 
                  DG_FP *u_y_avg, DG_FP *u_z_avg, DG_FP *v_jump, DG_FP *v_x_avg, 
                  DG_FP *v_y_avg, DG_FP *v_z_avg, DG_FP *w_jump, DG_FP *w_x_avg, 
                  DG_FP *w_y_avg, DG_FP *w_z_avg) {
  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * u_jump[i * DG_NPF + j];
      const DG_FP _sum = nx[i] * u_x_avg[i * DG_NPF + j]
                       + ny[i] * u_y_avg[i * DG_NPF + j];
                       + nz[i] * u_z_avg[i * DG_NPF + j];
      u_jump[i * DG_NPF + j]  = sJ[i] * (tau[i] * _jump - _sum);
      u_x_avg[i * DG_NPF + j] = nx[i] * sJ[i] * -_jump;
      u_y_avg[i * DG_NPF + j] = ny[i] * sJ[i] * -_jump;
      u_z_avg[i * DG_NPF + j] = nz[i] * sJ[i] * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * v_jump[i * DG_NPF + j];
      const DG_FP _sum = nx[i] * v_x_avg[i * DG_NPF + j]
                       + ny[i] * v_y_avg[i * DG_NPF + j];
                       + nz[i] * v_z_avg[i * DG_NPF + j];
      v_jump[i * DG_NPF + j]  = sJ[i] * (tau[i] * _jump - _sum);
      v_x_avg[i * DG_NPF + j] = nx[i] * sJ[i] * -_jump;
      v_y_avg[i * DG_NPF + j] = ny[i] * sJ[i] * -_jump;
      v_z_avg[i * DG_NPF + j] = nz[i] * sJ[i] * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * w_jump[i * DG_NPF + j];
      const DG_FP _sum = nx[i] * w_x_avg[i * DG_NPF + j]
                       + ny[i] * w_y_avg[i * DG_NPF + j];
                       + nz[i] * w_z_avg[i * DG_NPF + j];
      w_jump[i * DG_NPF + j]  = sJ[i] * (tau[i] * _jump - _sum);
      w_x_avg[i * DG_NPF + j] = nx[i] * sJ[i] * -_jump;
      w_y_avg[i * DG_NPF + j] = ny[i] * sJ[i] * -_jump;
      w_z_avg[i * DG_NPF + j] = nz[i] * sJ[i] * -_jump;
    }
  }
}