inline void ins_3d_proj_cg(const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                           const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                           const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                           const DG_FP *J, const DG_FP *pen, const DG_FP *rhs0,
                           const DG_FP *rhs1, const DG_FP *rhs2, DG_FP *u,
                           DG_FP *v, DG_FP *w, int *cell_g, int *conv_g,
                           DG_FP *iter_g) {
  const DG_FP *mass_mat = &dg_Mass_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *dr_mat = &dg_Dr_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &dg_Ds_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dg_Dt_kernel[(DG_ORDER - 1) * DG_NP * DG_NP];

  DG_FP r0[DG_NP], r1[DG_NP], r2[DG_NP];
  DG_FP r0_old[DG_NP], r1_old[DG_NP], r2_old[DG_NP];
  DG_FP x0[DG_NP], x1[DG_NP], x2[DG_NP];
  DG_FP p0[DG_NP], p1[DG_NP], p2[DG_NP];
  DG_FP tmp0[DG_NP], tmp1[DG_NP], tmp2[DG_NP];
  DG_FP tmpDiv[DG_NP], tmpMass[DG_NP];
  const DG_FP max_iter = 100;
  const DG_FP tol = 1e-24;
  DG_FP iter = 0;
  DG_FP residual;

  // r_{0} = b - Ax_{0}
  // Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, out0, out1, tmp0, tmp1);
  for(int i = 0; i < DG_NP; i++) {
    tmp0[i] = rx[0] * u[i] + ry[0] * v[i] + rz[0] * w[i];
    tmp1[i] = sx[0] * u[i] + sy[0] * v[i] + sz[0] * w[i];
    tmp2[i] = tx[0] * u[i] + ty[0] * v[i] + tz[0] * w[i];
  }
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, dr_mat, DG_NP, tmp0, 0.0, tmpDiv);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, ds_mat, DG_NP, tmp1, 1.0, tmpDiv);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, dt_mat, DG_NP, tmp2, 1.0, tmpDiv);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, tmpDiv, 0.0, tmpMass);
  op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, dr_mat, DG_NP, tmpMass, 0.0, tmp0);
  op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, ds_mat, DG_NP, tmpMass, 0.0, tmp1);
  op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, dt_mat, DG_NP, tmpMass, 0.0, tmp2);
  for(int i = 0; i < DG_NP; i++) {
    DG_FP r = tmp0[i]; DG_FP s = tmp1[i]; DG_FP t = tmp2[i];
    tmp0[i] = rx[0] * r + sx[0] * s + tx[0] * t;
    tmp1[i] = ry[0] * r + sy[0] * s + ty[0] * t;
    tmp2[i] = rz[0] * r + sz[0] * s + tz[0] * t;
  }
  op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, u, *pen, tmp0);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, v, *pen, tmp1);
  op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, w, *pen, tmp2);

  residual = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    r0[i] = rhs0[i] - tmp0[i];
    r1[i] = rhs1[i] - tmp1[i];
    r2[i] = rhs2[i] - tmp2[i];
    p0[i] = r0[i];
    p1[i] = r1[i];
    p2[i] = r2[i];
    residual += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
  }

  while(residual > tol && iter < max_iter) {
    // alpha = r^T r / (p^T A p)
    // Ax(op_0, op_1, op_2, op_3, Mass, J, *pen, p0, p1, tmp0, tmp1);
    for(int i = 0; i < DG_NP; i++) {
      tmp0[i] = rx[0] * p0[i] + ry[0] * p1[i] + rz[0] * p2[i];
      tmp1[i] = sx[0] * p0[i] + sy[0] * p1[i] + sz[0] * p2[i];
      tmp2[i] = tx[0] * p0[i] + ty[0] * p1[i] + tz[0] * p2[i];
    }
    op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, dr_mat, DG_NP, tmp0, 0.0, tmpDiv);
    op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, ds_mat, DG_NP, tmp1, 1.0, tmpDiv);
    op2_in_kernel_gemv(false, DG_NP, DG_NP, 1.0, dt_mat, DG_NP, tmp2, 1.0, tmpDiv);
    op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, tmpDiv, 0.0, tmpMass);
    op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, dr_mat, DG_NP, tmpMass, 0.0, tmp0);
    op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, ds_mat, DG_NP, tmpMass, 0.0, tmp1);
    op2_in_kernel_gemv(true, DG_NP, DG_NP, 1.0, dt_mat, DG_NP, tmpMass, 0.0, tmp2);
    for(int i = 0; i < DG_NP; i++) {
      DG_FP r = tmp0[i]; DG_FP s = tmp1[i]; DG_FP t = tmp2[i];
      tmp0[i] = rx[0] * r + sx[0] * s + tx[0] * t;
      tmp1[i] = ry[0] * r + sy[0] * s + ty[0] * t;
      tmp2[i] = rz[0] * r + sz[0] * s + tz[0] * t;
    }
    op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, p0, *pen, tmp0);
    op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, p1, *pen, tmp1);
    op2_in_kernel_gemv(false, DG_NP, DG_NP, J[0], mass_mat, DG_NP, p2, *pen, tmp2);

    DG_FP tmp_alpha_0 = 0.0;
    DG_FP tmp_alpha_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_alpha_0 += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
      tmp_alpha_1 += p0[i] * tmp0[i] + p1[i] * tmp1[i] + p2[i] * tmp2[i];
    }
    DG_FP alpha = tmp_alpha_0 / tmp_alpha_1;

    // x_{k+1} = x_{k} + alpha p_{k}
    // r_{k+1} = r_{k} - alpha Ap_{k}
    residual = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      u[i] = u[i] + alpha * p0[i];
      v[i] = v[i] + alpha * p1[i];
      w[i] = w[i] + alpha * p2[i];
      r0_old[i] = r0[i];
      r1_old[i] = r1[i];
      r2_old[i] = r2[i];
      r0[i] = r0[i] - alpha * tmp0[i];
      r1[i] = r1[i] - alpha * tmp1[i];
      r2[i] = r2[i] - alpha * tmp2[i];
      residual += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
    }

    if(residual < tol)
      break;

    // beta = r_{k+1}^{T} r_{k+1} / r_{k}^{T} r_{k}
    DG_FP tmp_beta_0 = 0.0;
    DG_FP tmp_beta_1 = 0.0;
    for(int i = 0; i < DG_NP; i++) {
      tmp_beta_0 += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
      tmp_beta_1 += r0_old[i] * r0_old[i] + r1_old[i] * r1_old[i] + r2_old[i] * r2_old[i];
    }
    DG_FP beta = tmp_beta_0 / tmp_beta_1;

    // p_{k+1} = r_{k+1} + beta p_{k}
    for(int i = 0; i < DG_NP; i++) {
      p0[i] = r0[i] + beta * p0[i];
      p1[i] = r1[i] + beta * p1[i];
      p2[i] = r2[i] + beta * p2[i];
    }

    iter++;
  }

  *cell_g += 1;
  *iter_g += (DG_FP) iter;
  if(iter < max_iter && !isnan(residual)) {
    *conv_g += 1;
  }
}
