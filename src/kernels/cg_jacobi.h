inline void cg_jacobi(const DG_FP *u_diag, const DG_FP *v_diag, const DG_FP *u_in, 
                      const DG_FP *v_in, DG_FP *u, DG_FP *v) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] = u_in[i] / u_diag[i];
    v[i] = v_in[i] / v_diag[i];
  }
}