inline void update_Q(const double *dt, double *s, const double *rk0,
                     const double *rk1, const double *rk2) {
  for(int i = 0; i < 15; i++) {
    double add = (*dt) * (rk0[i]/ 6.0 + rk1[i] / 6.0 + 2.0 * rk2[i] / 3.0);
    s[i] = s[i] + add;
  }
}
