inline void ls_group_modal(const double *modal, double *q) {
  // Group modal coefficients using quadratic mean
  #if DG_ORDER == 4
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[5] * modal[5];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[6] * modal[6] + modal[9] * modal[9];
  q[2] = sqrt(q[2] / 3.0);
  q[3] = modal[3] * modal[3] + modal[7] * modal[7] + modal[10] * modal[10]
         + modal[12] * modal[12];
  q[3] = sqrt(q[3] / 4.0);
  q[4] = modal[4] * modal[4] + modal[8] * modal[8] + modal[11] * modal[11]
         + modal[13] * modal[13] + modal[14] * modal[14];
  q[4] = sqrt(q[4] / 5.0);
  #elif DG_ORDER == 3
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[4] * modal[4];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[5] * modal[5] + modal[7] * modal[7];
  q[2] = sqrt(q[2] / 3.0);
  q[3] = modal[3] * modal[3] + modal[6] * modal[6] + modal[8] * modal[8]
         + modal[9] * modal[9];
  q[3] = sqrt(q[3] / 4.0);
  #elif DG_ORDER == 2
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[3] * modal[3];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[4] * modal[4] + modal[5] * modal[5];
  q[2] = sqrt(q[2] / 3.0);
  #else
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[2] * modal[2];
  q[1] = sqrt(q[1] / 2.0);
  #endif

  // Skyline pessimization
  #if DG_ORDER == 4
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]); q[3] = fabs(q[3]);
  q[4] = fabs(q[4]);
  q[4] = fmax(q[3], q[4]);
  q[3] = fmax(q[3], q[4]);
  q[2] = fmax(q[2], q[3]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #elif DG_ORDER == 3
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]); q[3] = fabs(q[3]);
  q[3] = fmax(q[2], q[3]);
  q[2] = fmax(q[2], q[3]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #elif DG_ORDER == 2
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]);
  q[2] = fmax(q[1], q[2]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #else
  q[0] = fabs(q[0]); q[1] = fabs(q[1]);
  q[1] = fmax(q[0], q[1]);
  q[0] = fmax(q[0], q[1]);
  #endif
}
