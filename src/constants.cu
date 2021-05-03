#include "constants.h"

Constants::Constants() {
  // Cubature constants
  cudaMalloc((void**)&cubDr_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubDr_d, cubDr_g, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cubDs_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubDs_d, cubDs_g, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cubV_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubV_d, cubV_g, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cubVDr_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubVDr_d, cubVDr_g, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cubVDs_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubVDs_d, cubVDs_g, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cubW_d, 46 * sizeof(double));
  cudaMemcpy(cubW_d, cubW_d, 46 * sizeof(double), cudaMemcpyHostToDevice);
  // Grad constants
  cudaMalloc((void**)&Dr_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dr_d, Dr_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&Drw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Drw_d, Drw_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&Ds_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Ds_d, Ds_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&Dsw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dsw_d, Dsw_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  // Gauss constants
  cudaMalloc((void**)&gaussW_d, 7 * sizeof(double));
  cudaMemcpy(gaussW_d, gaussW_g, 7 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF0Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF0Dr_d, gF0Dr_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF0DrR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF0DrR_d, gF0DrR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF0Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF0Ds_d, gF0Ds_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF0DsR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF0DsR_d, gF0DsR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF1Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF1Dr_d, gF1Dr_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF1DrR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF1DrR_d, gF1DrR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF1Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF1Ds_d, gF1Ds_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF1DsR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF1DsR_d, gF1DsR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF2Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF2Dr_d, gF2Dr_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF2DrR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF2DrR_d, gF2DrR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF2Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF2Ds_d, gF2Ds_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gF2DsR_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gF2DsR_d, gF2DsR_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp0_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp0_d, gFInterp0_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp0R_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp0R_d, gFInterp0R_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp1_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp1_d, gFInterp1_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp1R_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp1R_d, gFInterp1R_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp2_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp2_d, gFInterp2_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gFInterp2R_d, 7 * 15 * sizeof(double));
  cudaMemcpy(gFInterp2R_d, gFInterp2R_g, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&gInterp_d, 21 * 15 * sizeof(double));
  cudaMemcpy(gInterp_d, gInterp_g, 21 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  // Other constants
  cudaMalloc((void**)&invMass_d, 15 * 15 * sizeof(double));
  cudaMemcpy(invMass_d, invMass_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&r_d, 15 * sizeof(double));
  cudaMemcpy(r_d, r_g, 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&s_d, 15 * sizeof(double));
  cudaMemcpy(s_d, s_g, 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&ones_d, 15 * sizeof(double));
  cudaMemcpy(ones_d, ones_g, 15 * sizeof(double), cudaMemcpyHostToDevice);

  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
}

Constants::~Constants() {
  // Cubature constants
  cudaFree(cubDr_d);
  cudaFree(cubDs_d);
  cudaFree(cubV_d);
  cudaFree(cubVDr_d);
  cudaFree(cubVDs_d);
  cudaFree(cubW_d);
  // Grad constants
  cudaFree(Dr_d);
  cudaFree(Drw_d);
  cudaFree(Ds_d);
  cudaFree(Dsw_d);
  // Gauss constants
  cudaFree(gaussW_d);
  cudaFree(gF0Dr_d);
  cudaFree(gF0DrR_d);
  cudaFree(gF0Ds_d);
  cudaFree(gF0DsR_d);
  cudaFree(gF1Dr_d);
  cudaFree(gF1DrR_d);
  cudaFree(gF1Ds_d);
  cudaFree(gF1DsR_d);
  cudaFree(gF2Dr_d);
  cudaFree(gF2DrR_d);
  cudaFree(gF2Ds_d);
  cudaFree(gF2DsR_d);
  cudaFree(gFInterp0_d);
  cudaFree(gFInterp0R_d);
  cudaFree(gFInterp1_d);
  cudaFree(gFInterp1R_d);
  cudaFree(gFInterp2_d);
  cudaFree(gFInterp2R_d);
  cudaFree(gInterp_d);
  // Other constants
  cudaFree(invMass_d);
  cudaFree(LIFT_d);
  cudaFree(MASS_d);
  cudaFree(r_d);
  cudaFree(s_d);
  cudaFree(ones_d);

  cublasDestroy(handle);
}
