inline void calc_h_explicitly_3d(const DG_FP *x, const DG_FP *y, const DG_FP *z, DG_FP *h) {
  // Volume of the tet
  DG_FP tmp0[3], tmp1[3], tmpcross[3];
  tmp0[0] = x[1] - x[3];
  tmp0[1] = y[1] - y[3];
  tmp0[2] = z[1] - z[3];
  tmp1[0] = x[2] - x[3];
  tmp1[1] = y[2] - y[3];
  tmp1[2] = z[2] - z[3];
  tmpcross[0] = tmp0[1] * tmp1[2] - tmp0[2] * tmp1[1];
  tmpcross[1] = tmp0[2] * tmp1[0] - tmp0[0] * tmp1[2];
  tmpcross[2] = tmp0[0] * tmp1[1] - tmp0[1] * tmp1[0];
  tmp0[0] = x[0] - x[3];
  tmp0[1] = y[0] - y[3];
  tmp0[2] = z[0] - z[3];
  const DG_FP vol = fabs(tmp0[0] * tmpcross[0] + tmp0[1] * tmpcross[1] + tmp0[2] * tmpcross[2]) / 6.0;
  // Edge lengths from node 0
  const DG_FP a = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) + (z[0] - z[1]) * (z[0] - z[1]));
  const DG_FP b = sqrt((x[0] - x[2]) * (x[0] - x[2]) + (y[0] - y[2]) * (y[0] - y[2]) + (z[0] - z[2]) * (z[0] - z[2]));
  const DG_FP c = sqrt((x[0] - x[3]) * (x[0] - x[3]) + (y[0] - y[3]) * (y[0] - y[3]) + (z[0] - z[3]) * (z[0] - z[3]));
  // Opposite edge lengths
  const DG_FP a1 = sqrt((x[2] - x[3]) * (x[2] - x[3]) + (y[2] - y[3]) * (y[2] - y[3]) + (z[2] - z[3]) * (z[2] - z[3]));
  const DG_FP b1 = sqrt((x[1] - x[3]) * (x[1] - x[3]) + (y[1] - y[3]) * (y[1] - y[3]) + (z[1] - z[3]) * (z[1] - z[3]));
  const DG_FP c1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]) + (z[1] - z[2]) * (z[1] - z[2]));

  const DG_FP p = 0.5 * (a * a1 + b * b1 + c * c1);

  *h = sqrt(p * (p - a * a1) * (p - b * b1) * (p - c * c1)) / (6.0 * vol);
}