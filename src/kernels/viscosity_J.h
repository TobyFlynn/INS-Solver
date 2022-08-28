inline void viscosity_J(const int *p, const double *J, const double *u, 
                        const double *v, double *u_tmp, double *v_tmp) {
  for(int i = 0; i < DG_NP; i++) {
    u_tmp[i] = u[i] * J[i];
    v_tmp[i] = v[i] * J[i];
  }
}