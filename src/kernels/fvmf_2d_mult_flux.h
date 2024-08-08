inline void fvmf_2d_mult_flux(const DG_FP *nx, const DG_FP *ny, const DG_FP *sJ, 
                              const DG_FP *tau, const DG_FP *factor, DG_FP *u_jump, 
                              DG_FP *u_x_avg, DG_FP *u_y_avg, DG_FP *v_jump, 
                              DG_FP *v_x_avg, DG_FP *v_y_avg) {
  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * u_jump[i * DG_NPF + j];
      const DG_FP _sum = nx[i] * u_x_avg[i * DG_NPF + j]
                       + ny[i] * u_y_avg[i * DG_NPF + j];
      u_jump[i * DG_NPF + j]  = sJ[i] * (tau[i] * _jump - _sum);
      const int factor_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + i * DG_NPF + j];
      const DG_FP fact = factor[factor_ind];
      u_x_avg[i * DG_NPF + j] = nx[i] * sJ[i] * fact * -_jump;
      u_y_avg[i * DG_NPF + j] = ny[i] * sJ[i] * fact * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * v_jump[i * DG_NPF + j];
      const DG_FP _sum = nx[i] * v_x_avg[i * DG_NPF + j]
                       + ny[i] * v_y_avg[i * DG_NPF + j];
      v_jump[i * DG_NPF + j]  = sJ[i] * (tau[i] * _jump - _sum);
      const int factor_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + i * DG_NPF + j];
      const DG_FP fact = factor[factor_ind];
      v_x_avg[i * DG_NPF + j] = nx[i] * sJ[i] * fact * -_jump;
      v_y_avg[i * DG_NPF + j] = ny[i] * sJ[i] * fact * -_jump;
    }
  }
}