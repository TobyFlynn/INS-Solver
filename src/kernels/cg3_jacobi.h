inline void cg3_jacobi(const DG_FP *u_diag, const DG_FP *v_diag, const DG_FP *w_diag, 
                       const DG_FP *u_in, const DG_FP *v_in, const DG_FP *w_in, 
                       DG_FP *u, DG_FP *v, DG_FP *w) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] = u_in[i] / u_diag[i];
    v[i] = v_in[i] / v_diag[i];
    w[i] = w_in[i] / w_diag[i];
  }
}