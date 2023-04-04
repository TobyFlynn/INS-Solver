inline void ins_3d_proj_setup_0(const int *faceNum, const DG_FP *fscale,
                                DG_FP **f0) {
  f0[0][faceNum[0]] += fscale[0];
  f0[1][faceNum[1]] += fscale[1];
}
