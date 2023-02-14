/*
inline void Ax(const DG_FP *op_0, const DG_FP *op_1, const DG_FP *op_2,
               const DG_FP *op_3, const DG_FP *Mass,
               const DG_FP *J, const DG_FP pen, const DG_FP *in0,
               const DG_FP *in1, DG_FP *out0, DG_FP *out1) {
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
*/

inline void project_2d_cg(const DG_FP *mass_, const DG_FP *J, const DG_FP *op_0,
                          const DG_FP *op_1, const DG_FP *op_2, const DG_FP *op_3,
                          const DG_FP *pen, const DG_FP *rhs0,
                          const DG_FP *rhs1, DG_FP *out0, DG_FP *out1,
                          int *cell_g, int *conv_g, DG_FP *iter_g) {
  // Get matrices
  const DG_FP *Mass = &mass_[(DG_ORDER - 1) * DG_NP * DG_NP];

  DG_FP r0[DG_NP], r1[DG_NP], r0_old[DG_NP], r1_old[DG_NP], x0[DG_NP], x1[DG_NP];
  DG_FP p0[DG_NP], p1[DG_NP], tmp0[DG_NP], tmp1[DG_NP];
  const DG_FP max_iter = 100;
  const DG_FP tol = 1e-24;
  DG_FP iter = 0;
  DG_FP residual;

  for(int i = 0; i < DG_NP; i++) {
    out0[i] = 0.0;
    out1[i] = 0.0;
  }

  // r_{0} = b - Ax_{0}
  // Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, out0, out1, tmp0, tmp1);
  for(int i = 0; i < DG_NP; i++) {
    tmp0[i] = 0.0;
    tmp1[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      tmp0[i] += Mass[ind] * out0[j] + *pen * (op_0[ind] * out0[j] + op_3[ind] * out1[j]);
      tmp1[i] += Mass[ind] * out1[j] + *pen * (op_2[ind] * out0[j] + op_1[ind] * out1[j]);
    }
    tmp0[i] *= J[i];
    tmp1[i] *= J[i];
  }
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
    // Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, p0, p1, tmp0, tmp1);
    for(int i = 0; i < DG_NP; i++) {
      tmp0[i] = 0.0;
      tmp1[i] = 0.0;
      for(int j = 0; j < DG_NP; j++) {
        int ind = i + j * DG_NP;
        tmp0[i] += Mass[ind] * p0[j] + *pen * (op_0[ind] * p0[j] + op_3[ind] * p1[j]);
        tmp1[i] += Mass[ind] * p1[j] + *pen * (op_2[ind] * p0[j] + op_1[ind] * p1[j]);
      }
      tmp0[i] *= J[i];
      tmp1[i] *= J[i];
    }
    DG_FP tmp_alpha_0 = 0.0;
    DG_FP tmp_alpha_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_alpha_0 += r0[i] * r0[i] + r1[i] * r1[i];
      tmp_alpha_1 += p0[i] * tmp0[i] + p1[i] * tmp1[i];
    }
    DG_FP alpha = tmp_alpha_0 / tmp_alpha_1;

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
    DG_FP tmp_beta_0 = 0.0;
    DG_FP tmp_beta_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_beta_0 += r0[i] * r0[i] + r1[i] * r1[i];
      tmp_beta_1 += r0_old[i] * r0_old[i] + r1_old[i] * r1_old[i];
    }
    DG_FP beta = tmp_beta_0 / tmp_beta_1;

    // p_{k+1} = r_{k+1} + beta p_{k}
    for(int i = 0; i < DG_NP; i++) {
      p0[i] = r0[i] + beta * p0[i];
      p1[i] = r1[i] + beta * p1[i];
    }

    iter++;
  }

  *cell_g += 1;
  *iter_g += (DG_FP) iter;
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
