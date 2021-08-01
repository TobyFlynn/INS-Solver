inline void save_values(const double *v_vals, double *c_vals) {
  #if DG_ORDER == 4
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[5]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[5] + v_vals[6]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[6]) / 3.0;
  c_vals[3]  = (v_vals[2] + v_vals[6] + v_vals[7]) / 3.0;
  c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[7]) / 3.0;
  c_vals[5]  = (v_vals[3] + v_vals[7] + v_vals[8]) / 3.0;
  c_vals[6]  = (v_vals[3] + v_vals[4] + v_vals[8]) / 3.0;
  c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[9]) / 3.0;
  c_vals[8]  = (v_vals[6] + v_vals[9] + v_vals[10]) / 3.0;
  c_vals[9]  = (v_vals[6] + v_vals[7] + v_vals[10]) / 3.0;
  c_vals[10] = (v_vals[7] + v_vals[10] + v_vals[11]) / 3.0;
  c_vals[11] = (v_vals[7] + v_vals[8] + v_vals[11]) / 3.0;
  c_vals[12] = (v_vals[9] + v_vals[10] + v_vals[12]) / 3.0;
  c_vals[13] = (v_vals[10] + v_vals[12] + v_vals[13]) / 3.0;
  c_vals[14] = (v_vals[10] + v_vals[11] + v_vals[13]) / 3.0;
  c_vals[15] = (v_vals[12] + v_vals[13] + v_vals[14]) / 3.0;
  #elif DG_ORDER == 3
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[4]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[4] + v_vals[5]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[5]) / 3.0;
  c_vals[3]  = (v_vals[2] + v_vals[5] + v_vals[6]) / 3.0;
  c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[6]) / 3.0;
  c_vals[5]  = (v_vals[4] + v_vals[5] + v_vals[7]) / 3.0;
  c_vals[6]  = (v_vals[5] + v_vals[7] + v_vals[8]) / 3.0;
  c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[8]) / 3.0;
  c_vals[8]  = (v_vals[7] + v_vals[8] + v_vals[9]) / 3.0;
  #elif DG_ORDER == 2
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[3]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[3] + v_vals[4]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[4]) / 3.0;
  c_vals[3]  = (v_vals[3] + v_vals[4] + v_vals[5]) / 3.0;
  #else
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[2]) / 3.0;
  #endif
}
