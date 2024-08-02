inline void vmf_2d_mult_flux(const int *order, const DG_FP *nx, const DG_FP *ny,
                             const DG_FP *sJ, const DG_FP *tau, DG_FP *u_jump, 
                             DG_FP *u_x_avg, DG_FP *u_y_avg, DG_FP *v_jump,
                             DG_FP *v_x_avg, DG_FP *v_y_avg) {
  const int dg_npf = DG_CONSTANTS[(*order - 1) * DG_NUM_CONSTANTS + 1];

  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < dg_npf; j++) {
      const DG_FP _jump = 0.5 * u_jump[i * dg_npf + j];
      const DG_FP _sum = nx[i] * u_x_avg[i * dg_npf + j]
                       + ny[i] * u_y_avg[i * dg_npf + j];
      u_jump[i * dg_npf + j]  = sJ[i] * (tau[i] * _jump - _sum);
      u_x_avg[i * dg_npf + j] = nx[i] * sJ[i] * -_jump;
      u_y_avg[i * dg_npf + j] = ny[i] * sJ[i] * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    for(int j = 0; j < dg_npf; j++) {
      const DG_FP _jump = 0.5 * v_jump[i * dg_npf + j];
      const DG_FP _sum = nx[i] * v_x_avg[i * dg_npf + j]
                       + ny[i] * v_y_avg[i * dg_npf + j];
      v_jump[i * dg_npf + j]  = sJ[i] * (tau[i] * _jump - _sum);
      v_x_avg[i * dg_npf + j] = nx[i] * sJ[i] * -_jump;
      v_y_avg[i * dg_npf + j] = ny[i] * sJ[i] * -_jump;
    }
  }
}