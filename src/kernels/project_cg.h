inline void Ax(const double *op_0, const double *op_1, const double *op_2,
               const double *op_3, const double *Mass,
               const double *J, const double pen, const double *in0,
               const double *in1, double *out0, double *out1) {
  // Add mass term and div-div term to output
  for(int i = 0; i < DG_NP; i++) {
    out0[i] = 0.0;
    out1[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      out0[i] += Mass[ind] * in0[j] + pen * (op_0[ind] * in0[j] + op_3[ind] * in1[j]);
      out1[i] += Mass[ind] * in1[j] + pen * (op_2[ind] * in0[j] + op_1[ind] * in1[j]);
      // out0[i] += Mass[ind] * in0[j] - pen * (op_0[ind] * in0[j] + op_2[ind] * in1[j]);
      // out1[i] += Mass[ind] * in1[j] - pen * (op_3[ind] * in0[j] + op_1[ind] * in1[j]);

      // out0[i] += Mass[ind] * in0[j];
      // out1[i] += Mass[ind] * in1[j];
    }
    out0[i] *= J[i];
    out1[i] *= J[i];
  }
}

inline void project_cg(const double *mass_, const double *J, const double *op_0,
                       const double *op_1, const double *op_2, const double *op_3,
                       const double *pen, const double *rhs0,
                       const double *rhs1, double *out0, double *out1,
                       int *cell_g, int *conv_g, double *iter_g) {
  // Get matrices
  const double *Mass = &mass_[(DG_ORDER - 1) * DG_NP * DG_NP];

  double r0[DG_NP], r1[DG_NP], r0_old[DG_NP], r1_old[DG_NP], x0[DG_NP], x1[DG_NP];
  double p0[DG_NP], p1[DG_NP], tmp0[DG_NP], tmp1[DG_NP];
  const double max_iter = 100;
  const double tol = 1e-24;
  double iter = 0;
  double residual;

  for(int i = 0; i < DG_NP; i++) {
    out0[i] = 0.0;
    out1[i] = 0.0;
  }

  // r_{0} = b - Ax_{0}
  Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, out0, out1, tmp0, tmp1);
  residual = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    r0[i] = rhs0[i] - tmp0[i];
    r1[i] = rhs1[i] - tmp1[i];
    p0[i] = r0[i];
    p1[i] = r1[i];
    residual += r0[i] * r0[i] + r1[i] * r1[i];
  }

  while(residual > tol && iter < max_iter) {
    // alpha = r^T r / (p^T A p)
    Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, p0, p1, tmp0, tmp1);
    double tmp_alpha_0 = 0.0;
    double tmp_alpha_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_alpha_0 += r0[i] * r0[i] + r1[i] * r1[i];
      tmp_alpha_1 += p0[i] * tmp0[i] + p1[i] * tmp1[i];
    }
    double alpha = tmp_alpha_0 / tmp_alpha_1;

    // x_{k+1} = x_{k} + alpha p_{k}
    // r_{k+1} = r_{k} - alpha Ap_{k}
    residual = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      out0[i] = out0[i] + alpha * p0[i];
      out1[i] = out1[i] + alpha * p1[i];
      r0_old[i] = r0[i];
      r1_old[i] = r1[i];
      r0[i] = r0[i] - alpha * tmp0[i];
      r1[i] = r1[i] - alpha * tmp1[i];
      residual += r0[i] * r0[i] + r1[i] * r1[i];
    }

    if(residual < tol)
      break;

    // beta = r_{k+1}^{T} r_{k+1} / r_{k}^{T} r_{k}
    double tmp_beta_0 = 0.0;
    double tmp_beta_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_beta_0 += r0[i] * r0[i] + r1[i] * r1[i];
      tmp_beta_1 += r0_old[i] * r0_old[i] + r1_old[i] * r1_old[i];
    }
    double beta = tmp_beta_0 / tmp_beta_1;

    // p_{k+1} = r_{k+1} + beta p_{k}
    for(int i = 0; i < DG_NP; i++) {
      p0[i] = r0[i] + beta * p0[i];
      p1[i] = r1[i] + beta * p1[i];
    }

    iter++;
  }

  *cell_g += 1;
  *iter_g += (double) iter;
  if(iter < max_iter && !isnan(residual)) {
    *conv_g += 1;
  } else {
    printf("%g\n", residual);
    for(int i = 0; i < DG_NP; i++) {
      out0[i] = -2.0;
      out1[i] = -2.0;
    }
  }
}
