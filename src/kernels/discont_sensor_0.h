inline void discont_sensor_0(const DG_FP *simplexes, const DG_FP *in,
                             const DG_FP *u_modal, DG_FP *u_hat) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP *simplex_ptr = simplexes + i * DG_NP;
    DG_FP tmp = 0.0;
    for(int j = 0; j < DG_NP; j++) {
      tmp += simplex_ptr[j] * u_modal[j];
    }
    u_hat[i] = in[i] - tmp;
  }
}
