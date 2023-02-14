inline void calc_h_ls(const DG_FP *x, const DG_FP *y, const DG_FP *s, DG_FP *dt) {
  DG_FP len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  DG_FP len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  DG_FP len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));
  DG_FP sper = (len0 + len1 + len2) / 2.0;
  DG_FP area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  if(*dt < area / sper && fabs(s[0]) < 0.04)
    *dt = area / sper;
}
