inline void ins_3d_advec_0(const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *f00, DG_FP *f01, DG_FP *f02, DG_FP *f10,
                           DG_FP *f11, DG_FP *f12, DG_FP *f20, DG_FP *f21,
                           DG_FP *f22) {
  for(int i = 0; i < DG_NP; i++) {
    const DG_FP _u = u[i];
    const DG_FP _v = v[i];
    const DG_FP _w = w[i];
    f00[i] = _u * _u; f01[i] = _u * _v; f02[i] = _u * _w;
    f10[i] = _v * _u; f11[i] = _v * _v; f12[i] = _v * _w;
    f20[i] = _w * _u; f21[i] = _w * _v; f22[i] = _w * _w;
  }
}
