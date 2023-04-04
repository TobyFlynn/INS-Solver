inline void pmf_3d_mult_faces_flux(const int *p, const int *faceNums,
                          const int *fmaskF, const DG_FP *nx, const DG_FP *ny,
                          const DG_FP *nz, const DG_FP *fscale, const DG_FP *sJ,
                          const DG_FP *in, const DG_FP **in_p, const DG_FP *in_x,
                          const DG_FP *in_y, const DG_FP *in_z, const DG_FP **in_x_p,
                          const DG_FP **in_y_p, const DG_FP **in_z_p, DG_FP *l_x,
                          DG_FP *l_y, DG_FP *l_z, DG_FP *out) {
  const int dg_np  = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS];
  const int dg_npf = DG_CONSTANTS[(*p - 1) * DG_NUM_CONSTANTS + 1];
  const int *fmask  = &FMASK[(*p - 1) * 4 * DG_NPF];

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[i * 2], fscale[i * 2 + 1]);
    const DG_FP int_fact = 0.5 * sJ[i];
    const DG_FP *inR = in_p[i];
    const DG_FP *in_xR = in_x_p[i];
    const DG_FP *in_yR = in_y_p[i];
    const DG_FP *in_zR = in_z_p[i];

    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];

    #pragma omp simd
    for(int j = 0; j < dg_npf; j++) {
      const int fmaskIndL = fmaskL[j];
      const int fmaskIndR = fmaskR[j];
      const DG_FP diffL_u = in[fmaskIndL] - inR[fmaskIndR];
      const DG_FP diffL_u_x = nx[i] * (in_xR[fmaskIndR] + in_x[fmaskIndL]);
      const DG_FP diffL_u_y = ny[i] * (in_yR[fmaskIndR] + in_y[fmaskIndL]);
      const DG_FP diffL_u_z = nz[i] * (in_zR[fmaskIndR] + in_z[fmaskIndL]);
      const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

      const int indL = findL + j;
      out[indL] = int_fact * (gtau * diffL_u - diffL_u_grad);
      const DG_FP l_tmpL = int_fact * -diffL_u;
      l_x[indL] = nx[i] * l_tmpL;
      l_y[indL] = ny[i] * l_tmpL;
      l_z[indL] = nz[i] * l_tmpL;
    }
  }
}
