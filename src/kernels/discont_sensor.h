inline void discont_sensor(const DG_FP *e0, const DG_FP *s0, const DG_FP *k,
                           const DG_FP *geof, const DG_FP *u,
                           const DG_FP *u_hat, DG_FP *out) {
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
  const DG_FP c = 1.0;
  DG_FP se = log10(fmin(c * DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER * (u_hat_ip / u_ip), 1.0));
  // printf("%g %g - s0 = %g   k = %g\n", u_hat_ip / u_ip, se, *s0, *k);
  if(se < *s0 - *k) {
    *out = 0.0;
  } else if(se > *s0 + *k) {
    *out = *e0;
  } else {
    *out = (*e0 / 2.0) * (1.0 + sin(PI * (se - *s0) / (*k * 2.0)));
  }
}
