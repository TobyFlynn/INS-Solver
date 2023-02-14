inline void discont_sensor(const DG_FP *e0, const DG_FP *s0, const DG_FP *k,
                           const DG_FP *mass, const DG_FP *J, const DG_FP *u,
                           const DG_FP *u_hat, DG_FP *out) {
  const DG_FP *mat = &mass[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP PI = 3.141592653589793238463;
  DG_FP tmp_u[DG_NP];
  DG_FP tmp_u_hat[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    tmp_u[i] = 0.0;
    tmp_u_hat[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      tmp_u[i] += mat[ind] * u[j];
      tmp_u_hat[i] += mat[ind] * u_hat[j];
    }
  }

  DG_FP u_ip = 0.0;
  DG_FP u_hat_ip = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    u_ip += *J * tmp_u[i] * u[i];
    u_hat_ip += *J * tmp_u_hat[i] * u_hat[i];
  }

  DG_FP se = log10(u_hat_ip / u_ip);
  // printf("%g %g - s0 = %g   k = %g\n", u_hat_ip / u_ip, se, *s0, *k);
  if(se < *s0 - *k) {
    *out = 0.0;
  } else if(se > *s0 + *k) {
    *out = *e0;
  } else {
    *out = (*e0 / 2.0) * (1.0 + sin(PI * (se - *s0) / (*k * 2.0)));
  }
}