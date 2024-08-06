inline void fvmf_3d_mult_flux(const DG_FP *nx, const DG_FP *ny, const DG_FP *nz, 
                              const DG_FP *sJ, const DG_FP *tau, const DG_FP *factor, 
                              DG_FP *u_jump, DG_FP *u_avg_x, DG_FP *u_avg_y, DG_FP *u_avg_z,
                              DG_FP *v_jump, DG_FP *v_avg_x, DG_FP *v_avg_y, DG_FP *v_avg_z,
                              DG_FP *w_jump, DG_FP *w_avg_x, DG_FP *w_avg_y, DG_FP *w_avg_z) {
  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP _nx = nx[i];
    const DG_FP _ny = ny[i];
    const DG_FP _nz = nz[i];
    const DG_FP _sJ = sJ[i];
    const DG_FP _tau = tau[i];
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * u_jump[i * DG_NPF + j];
      const DG_FP _sum = _nx * u_avg_x[i * DG_NPF + j]
                       + _ny * u_avg_y[i * DG_NPF + j]
                       + _nz * u_avg_z[i * DG_NPF + j];
      u_jump[i * DG_NPF + j]  = _sJ * (_tau * _jump - _sum);
      const int factor_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + i * DG_NPF + j];
      const DG_FP fact = factor[factor_ind];
      u_avg_x[i * DG_NPF + j] = _nx * _sJ * fact * -_jump;
      u_avg_y[i * DG_NPF + j] = _ny * _sJ * fact * -_jump;
      u_avg_z[i * DG_NPF + j] = _nz * _sJ * fact * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP _nx = nx[i];
    const DG_FP _ny = ny[i];
    const DG_FP _nz = nz[i];
    const DG_FP _sJ = sJ[i];
    const DG_FP _tau = tau[i];
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * v_jump[i * DG_NPF + j];
      const DG_FP _sum = _nx * v_avg_x[i * DG_NPF + j]
                       + _ny * v_avg_y[i * DG_NPF + j]
                       + _nz * v_avg_z[i * DG_NPF + j];
      v_jump[i * DG_NPF + j]  = _sJ * (_tau * _jump - _sum);
      const int factor_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + i * DG_NPF + j];
      const DG_FP fact = factor[factor_ind];
      v_avg_x[i * DG_NPF + j] = _nx * _sJ * fact * -_jump;
      v_avg_y[i * DG_NPF + j] = _ny * _sJ * fact * -_jump;
      v_avg_z[i * DG_NPF + j] = _nz * _sJ * fact * -_jump;
    }
  }

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP _nx = nx[i];
    const DG_FP _ny = ny[i];
    const DG_FP _nz = nz[i];
    const DG_FP _sJ = sJ[i];
    const DG_FP _tau = tau[i];
    for(int j = 0; j < DG_NPF; j++) {
      const DG_FP _jump = 0.5 * w_jump[i * DG_NPF + j];
      const DG_FP _sum = _nx * w_avg_x[i * DG_NPF + j]
                       + _ny * w_avg_y[i * DG_NPF + j]
                       + _nz * w_avg_z[i * DG_NPF + j];
      w_jump[i * DG_NPF + j]  = _sJ * (_tau * _jump - _sum);
      const int factor_ind = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + i * DG_NPF + j];
      const DG_FP fact = factor[factor_ind];
      w_avg_x[i * DG_NPF + j] = _nx * _sJ * fact * -_jump;
      w_avg_y[i * DG_NPF + j] = _ny * _sJ * fact * -_jump;
      w_avg_z[i * DG_NPF + j] = _nz * _sJ * fact * -_jump;
    }
  }
}