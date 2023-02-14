inline void euler_2d_ic(const DG_FP *x, const DG_FP *y, DG_FP *q0,
                        DG_FP *q1, DG_FP *q2, DG_FP *q3) {
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP R = 1.5;
  const DG_FP S = 13.5;
  const DG_FP M = 0.4;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP f = (1.0 - x[i] * x[i] - y[i] * y[i]) / (2.0 * R * R);
    DG_FP u = (S * y[i] * exp(f)) / (2.0 * PI * R);
    DG_FP v = (1.0 - ((S * x[i] * exp(f)) / (2.0 * PI * R)));
    DG_FP rho_tmp = 1.0 - ((S * S * M * M * (gamma_e - 1.0) * exp(2.0 * f)) / (8.0 * PI * PI));
    DG_FP rho = pow(rho_tmp, 1.0 / (gamma_e - 1.0));
    DG_FP p = pow(rho, gamma_e) / (gamma_e * M * M);

    q0[i] = rho;
    q1[i] = rho * u;
    q2[i] = rho * v;
    q3[i] = 0.5 * (rho * u * u + rho * v * v) + p / (gamma_e - 1.0);
  }
}