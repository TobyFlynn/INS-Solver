inline void fpmf_2d_mult_faces(const int **p, const int *edgeNum,
                  const bool *reverse, const DG_FP **sJ, const DG_FP **nx,
                  const DG_FP **ny, const DG_FP **h, const DG_FP **factor,
                  const DG_FP **in, const DG_FP **in_x, const DG_FP **in_y,
                  DG_FP **out, DG_FP **l_x, DG_FP **l_y) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  // Get constants
  // Using same Gauss points so should be able to replace dg_gf_npL and
  // dg_gf_npR with DG_GF_NP
  const int dg_npL      = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npfL     = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_npL   = DG_CONSTANTS[(p[0][0] - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussWL = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];
  const int dg_npR      = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS];
  const int dg_npfR     = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS + 1];
  const int dg_gf_npR   = DG_CONSTANTS[(p[1][0] - 1) * DG_NUM_CONSTANTS + 4];
  const DG_FP *gaussWR = &gaussW_g[(p[1][0] - 1) * DG_GF_NP];

  // Left edge
  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;

  // Left edge
  DG_FP tau[DG_GF_NP];
  DG_FP maxtau = 0.0;
  DG_FP max_hinv = fmax(h[0][edgeL * dg_npfL], h[1][edgeR * dg_npfR]);
  for(int i = 0; i < DG_GF_NP; i++) {
    int indL = exIndL + i;
    int indR;
    if(reverse)
      indR = exIndR + DG_GF_NP - 1 - i;
    else
      indR = exIndR + i;

    tau[i] = 0.5 * max_hinv * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(factor[0][indL], factor[1][indR]);
    if(tau[i] > maxtau) maxtau = tau[i];
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    tau[i] = maxtau;
  }

  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd = exIndL + i;
    int rInd;
    int rWInd;
    if(*reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
      rWInd = DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
      rWInd = i;
    }

    const DG_FP diffL_u = in[0][lInd] - in[1][rInd];
    const DG_FP diffL_u_x = nx[0][lInd] * (factor[1][rInd] * in_x[1][rInd] + factor[0][lInd] * in_x[0][lInd]);
    const DG_FP diffL_u_y = ny[0][lInd] * (factor[1][rInd] * in_y[1][rInd] + factor[0][lInd] * in_y[0][lInd]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y;

    out[0][lInd] += 0.5 * gaussWL[i] * sJ[0][lInd] * (tau[i] * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * factor[0][lInd] * gaussWL[i] * sJ[0][lInd] * -diffL_u;
    l_x[0][lInd] += nx[0][lInd] * l_tmpL;
    l_y[0][lInd] += ny[0][lInd] * l_tmpL;

    const DG_FP diffR_u = in[1][rInd] - in[0][lInd];
    const DG_FP diffR_u_x = nx[1][rInd] * (factor[1][rInd] * in_x[1][rInd] + factor[0][lInd] * in_x[0][lInd]);
    const DG_FP diffR_u_y = ny[1][rInd] * (factor[1][rInd] * in_y[1][rInd] + factor[0][lInd] * in_y[0][lInd]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y;

    out[1][rInd] += 0.5 * gaussWR[rWInd] * sJ[1][rInd] * (tau[rWInd] * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * factor[1][rInd] * gaussWR[rWInd] * sJ[1][rInd] * -diffR_u;
    l_x[1][rInd] += nx[1][rInd] * l_tmpR;
    l_y[1][rInd] += ny[1][rInd] * l_tmpR;
  }
}
