inline void calc_min_h_2d(const double *x, const double *y, double *h) {
  double len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  double len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  double len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));
  double sper = (len0 + len1 + len2) / 2.0;
  double area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  if(*h > area / sper)
    *h = area / sper;
}
