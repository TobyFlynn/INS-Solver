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
  cudaMalloc((void**)&lift_d, 15 * 15 * sizeof(double));
  cudaMemcpy(lift_d, lift_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&mass_d, 15 * 15 * sizeof(double));
  cudaMemcpy(mass_d, mass_g, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
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
  cudaFree(lift_d);
  cudaFree(mass_d);
  cudaFree(r_d);
  cudaFree(s_d);
  cudaFree(ones_d);

  cublasDestroy(handle);
}

double* Constants::get_ptr(Constant_Matrix mat) {
  switch(mat) {
    case CUB_DR:
      return cubDr_d;
    case CUB_DS:
      return cubDs_d;
    case CUB_V:
      return cubV_d;
    case CUB_VDR:
      return cubVDr_d;
    case CUB_VDS:
      return cubVDs_d;
    case CUB_W:
      return cubW_d;
    case DR:
      return Dr_d;
    case DRW:
      return Drw_d;
    case DS:
      return Ds_d;
    case DSW:
      return Dsw_d;
    case GAUSS_W:
      return gaussW_d;
    case GAUSS_F0DR:
      return gF0Dr_d;
    case GAUSS_F0DR_R:
      return gF0DrR_d;
    case GAUSS_F0DS:
      return gF0Ds_d;
    case GAUSS_F0DS_R:
      return gF0DsR_d;
    case GAUSS_F1DR:
      return gF1Dr_d;
    case GAUSS_F1DR_R:
      return gF1DrR_d;
    case GAUSS_F1DS:
      return gF1Ds_d;
    case GAUSS_F1DS_R:
      return gF1DsR_d;
    case GAUSS_F2DR:
      return gF2Dr_d;
    case GAUSS_F2DR_R:
      return gF2DrR_d;
    case GAUSS_F2DS:
      return gF2Ds_d;
    case GAUSS_F2DS_R:
      return gF2DsR_d;
    case GAUSS_FINTERP0:
      return gFInterp0_d;
    case GAUSS_FINTERP0_R:
      return gFInterp0R_d;
    case GAUSS_FINTERP1:
      return gFInterp1_d;
    case GAUSS_FINTERP1_R:
      return gFInterp1R_d;
    case GAUSS_FINTERP2:
      return gFInterp2_d;
    case GAUSS_FINTERP2_R:
      return gFInterp2R_d;
    case GAUSS_INTERP:
      return gInterp_d;
    case INV_MASS:
      return invMass_d;
    case LIFT:
      return lift_d;
    case MASS:
      return mass_d;
    case R:
      return r_d;
    case S:
      return s_d;
    case ONES:
      return ones_d;
  }
  return nullptr;
}
