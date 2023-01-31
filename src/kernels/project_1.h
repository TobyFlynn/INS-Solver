inline void project_1(const double *mass_, const double *J, const double *op_0,
                      const double *op_1, const double *pen, const double *in0,
                      const double *in1, double *out0, double *out1) {
  // Get matrices
  const double *Mass = &mass_[(DG_ORDER - 1) * DG_NP * DG_NP];

  // Add mass term and div-div term to output
  for(int i = 0; i < DG_NP; i++) {
    out0[i] = 0.0;
    out1[i] = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      int ind = i + j * DG_NP;
      out0[i] += Mass[ind] * in0[j] + (*pen) * (op_0[ind] * in0[j] + op_1[ind] * in1[j]);
      out1[i] += Mass[ind] * in1[j] + (*pen) * (op_0[ind] * in0[j] + op_1[ind] * in1[j]);
    }
    out0[i] *= J[i];
    out1[i] *= J[i];
  }
}
