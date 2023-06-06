inline void ins_proj_setup_0(const int *faceNum, const DG_FP *fscale, DG_FP **f0) {
  const int faceNumL = faceNum[0];
  const int faceNumR = faceNum[1];
  f0[0][faceNumL] += fscale[0];
  f0[1][faceNumR] += fscale[1];
}
