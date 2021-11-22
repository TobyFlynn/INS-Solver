inline void ls_rhs(const double *sign, const double *dpldx, const double *dprdx,
                   const double *dpldy, const double *dprdy, double *rk) {
  for(int i = 0; i < DG_NP; i++) {
    if(sign[i] > 0.0) {
      double plmx2 = fmin(dpldx[i], 0.0);
      plmx2 = plmx2 * plmx2;
      double plmy2 = fmin(dpldy[i], 0.0);
      plmy2 = plmy2 * plmy2;
      double prpx2 = fmax(dprdx[i], 0.0);
      prpx2 = prpx2 * prpx2;
      double prpy2 = fmax(dprdy[i], 0.0);
      prpy2 = prpy2 * prpy2;

      double H = sign[i] * (sqrt(fmax(plmx2, prpx2) + fmax(plmy2, prpy2)) - 1.0);
      rk[i] = -H;
    } else {
      double plpx2 = fmax(dpldx[i], 0.0);
      plpx2 = plpx2 * plpx2;
      double plpy2 = fmax(dpldy[i], 0.0);
      plpy2 = plpy2 * plpy2;
      double prmx2 = fmin(dprdx[i], 0.0);
      prmx2 = prmx2 * prmx2;
      double prmy2 = fmin(dprdy[i], 0.0);
      prmy2 = prmy2 * prmy2;

      double H = sign[i] * (sqrt(fmax(plpx2, prmx2) + fmax(plpy2, prmy2)) - 1.0);
      rk[i] = -H;
    }
  }
}
