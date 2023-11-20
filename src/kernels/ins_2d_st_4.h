inline void ins_2d_st_4(const DG_FP *alpha_, const DG_FP *nx, const DG_FP *ny, 
                        const DG_FP *sJ, const DG_FP *geof, DG_FP *mS, DG_FP *pS) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    const int face_ind = i/DG_CUB_SURF_2D_NP;
    const DG_FP _nx = nx[face_ind];
    const DG_FP _ny = ny[face_ind];
    const DG_FP _fscale = sJ[face_ind] / geof[J_IND];
    // const DG_FP _fscale = sJ[face_ind];
    // const DG_FP _mS = 0.5 * tanh(PI * mS[i] / alpha);
    // const DG_FP _pS = 0.5 * tanh(PI * pS[i] / alpha);
    DG_FP _mS = 0.5 * (1 + mS[i] / alpha + sin(PI * mS[i] / alpha) / PI);
    if(mS[i] < -alpha) _mS = 0.0;
    if(mS[i] > alpha) _mS = 1.0;
    DG_FP _pS = 0.5 * (1 + pS[i] / alpha + sin(PI * pS[i] / alpha) / PI);
    if(pS[i] < -alpha) _pS = 0.0;
    if(pS[i] > alpha) _pS = 1.0;

    mS[i] = _fscale * _nx * 0.5 * (_mS + _pS);
    pS[i] = _fscale * _ny * 0.5 * (_mS + _pS);
  }
}