// Simplistic kernel that works for the problem and BC here but not in general
inline void set_rhs(const double *J, const double *uD, double *rhs) {
  for(int i = 0; i < 15; i++) {
    rhs[FMASK[i]] += uD[i];
  }
  for(int i = 0; i < 15; i++) {
    rhs[i] *= J[i];
  }
}
