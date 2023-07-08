inline void discont_sensor_filter(const DG_FP *_max_alpha, const DG_FP *s0,
                                  const DG_FP *k, const DG_FP *c,
                                  const DG_FP *geof, const DG_FP *u,
                                  const DG_FP *u_hat, DG_FP *modal) {
  const DG_FP *mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP PI = 3.141592653589793238463;
  DG_FP tmp[DG_NP];
  DG_FP tmp_u[DG_NP];
  DG_FP tmp_u_hat[DG_NP];

  for(int i = 0; i < DG_NP; i++) {
    tmp[i] = u[i] * u[i] * geof[J_IND];
  }

  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, mat, DG_NP, tmp, 0.0, tmp_u);

  for(int i = 0; i < DG_NP; i++) {
    tmp[i] = u_hat[i] * u_hat[i] * geof[J_IND];
  }

  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, mat, DG_NP, tmp, 0.0, tmp_u_hat);

  DG_FP u_ip = 0.0;
  DG_FP u_hat_ip = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    u_ip += tmp_u[i];
    u_hat_ip += tmp_u_hat[i];
  }

  // DG_FP se = log10(u_hat_ip / u_ip);
  DG_FP se = log10(fmin(*c * DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER * (u_hat_ip / u_ip), 1.0));
  const DG_FP min_alpha = 0.0;
  const DG_FP max_alpha = *_max_alpha;
  const int cutoff_order = DG_ORDER;
  const int filter_strength = 32;

  if(se < *s0 - *k) {
    // No filtering
  } else if(se > *s0 + *k) {
    // Max filtering
    int modal_ind = 0;
    for(int i = 0; i < DG_ORDER + 1; i++) {
      for(int j = 0; j < DG_ORDER - i + 1; j++) {
        for(int k = 0; k < DG_ORDER - i - j + 1; k++) {
          if(i + j + k >= cutoff_order) {
            modal[modal_ind] *= exp(-max_alpha * pow(((i + j + k - cutoff_order + 1)/(DG_ORDER - cutoff_order + 1)),filter_strength));
          }
          modal_ind++;
        }
      }
    }
  } else {
    // Some filtering
    const DG_FP local_alpha = min_alpha + ((max_alpha - min_alpha) / 2.0) * (1.0 + sin(PI * (se - *s0) / (*k * 2.0)));
    int modal_ind = 0;
    for(int i = 0; i < DG_ORDER + 1; i++) {
      for(int j = 0; j < DG_ORDER - i + 1; j++) {
        for(int k = 0; k < DG_ORDER - i - j + 1; k++) {
          if(i + j + k >= cutoff_order) {
            modal[modal_ind] *= exp(-local_alpha * pow(((i + j + k - cutoff_order + 1)/(DG_ORDER - cutoff_order + 1)),filter_strength));
          }
          modal_ind++;
        }
      }
    }
  }
}
