inline void calc_h_explicitly(const DG_FP *x, const DG_FP *y, DG_FP *h) {
  DG_FP len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  DG_FP len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  DG_FP len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));

  DG_FP s = 0.5 * (len0 + len1 + len2);
  DG_FP radius = len0 * len1 * len2;
  radius = radius / (4.0 * sqrt(s * (s - len0) * (s - len1) * (s - len2)));
  *h = radius;
  // DG_FP sper = (len0 + len1 + len2) / 2.0;
  // DG_FP area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  // *h = area / sper;
}
