inline void ins_vis_bc_2d(const DG_FP *t, const int *bedge_type,
                          const int *bedgeNum, const DG_FP *x, const DG_FP *y,
                          const DG_FP *nx, const DG_FP *ny, DG_FP *exQ0,
                          DG_FP *exQ1) {
  int exInd = *bedgeNum * DG_GF_NP;
  const DG_FP PI = 3.141592653589793238463;

  if(*bedge_type == 0) {
    // Inflow - BC function dependant on time
    // for(int i = 0; i < DG_GF_NP; i++) {
    //   exQ0[exInd + i] += 1.0;
    // }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
  } else {
    // Wall - No slip
  }
}
