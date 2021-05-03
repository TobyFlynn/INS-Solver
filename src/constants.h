#ifndef __CONSTANTS_H
#define __CONSTANTS_H

extern double cubDr_g[46*15];
extern double cubDs_g[46*15];
extern double cubV_g[46*15];
extern double cubVDr_g[46*15];
extern double cubVDs_g[46*15];
extern double cubW_g[46];

extern double Dr_g[15*15];
extern double Drw_g[15*15];
extern double Ds_g[15*15];
extern double Dsw_g[15*15];

extern double gaussW_g[7];
extern double gF0Dr_g[7*15];
extern double gF0DrR_g[7*15];
extern double gF0Ds_g[7*15];
extern double gF0DsR_g[7*15];
extern double gF1Dr_g[7*15];
extern double gF1DrR_g[7*15];
extern double gF1Ds_g[7*15];
extern double gF1DsR_g[7*15];
extern double gF2Dr_g[7*15];
extern double gF2DrR_g[7*15];
extern double gF2Ds_g[7*15];
extern double gF2DsR_g[7*15];
extern double gFInterp0_g[7*15];
extern double gFInterp0R_g[7*15];
extern double gFInterp1_g[7*15];
extern double gFInterp1R_g[7*15];
extern double gFInterp2_g[7*15];
extern double gFInterp2R_g[7*15];
extern double gInterp_g[21*15];

extern double invMass_g[15*15];
extern double lift_g[15*15];
extern double mass_g[15*15];
extern double r_g[15];
extern double s_g[15];
extern double ones_g[15];

#ifdef INS_CUDA
#include "cublas_v2.h"
#endif

class Constants {
public:
  enum Constant_Matrix {
    CUB_DR, CUB_DS, CUB_V, CUB_VDR, CUB_VDS, CUB_W,

    DR, DRW, DS, DSW,

    GAUSS_W, GAUSS_F0DR, GAUSS_F0DR_R, GAUSS_F0DS, GAUSS_F0DS_R,
    GAUSS_F1DR, GAUSS_F1DR_R, GAUSS_F1DS, GAUSS_F1DS_R,
    GAUSS_F2DR, GAUSS_F2DR_R, GAUSS_F2DS, GAUSS_F2DS_R,
    GAUSS_FINTERP0, GAUSS_FINTERP0_R, GAUSS_FINTERP1, GAUSS_FINTERP1_R,
    GAUSS_FINTERP2, GAUSS_FINTERP2_R, GAUSS_INTERP,

    INV_MASS, LIFT, MASS, R, S, ONES
  };

  Constants();
  ~Constants();

  double* get_ptr(Constant_Matrix mat);

  double *cubDr, *cubDr_d;
  double *cubDs, *cubDs_d;
  double *cubV, *cubV_d;
  double *cubVDr, *cubVDr_d;
  double *cubVDs, *cubVDs_d;
  double *cubW, *cubW_d;

  double *Dr, *Dr_d;
  double *Drw, *Drw_d;
  double *Ds, *Ds_d;
  double *Dsw, *Dsw_d;

  double *gaussW, *gaussW_d;
  double *gF0Dr, *gF0Dr_d;
  double *gF0DrR, *gF0DrR_d;
  double *gF0Ds, *gF0Ds_d;
  double *gF0DsR, *gF0DsR_d;
  double *gF1Dr, *gF1Dr_d;
  double *gF1DrR, *gF1DrR_d;
  double *gF1Ds, *gF1Ds_d;
  double *gF1DsR, *gF1DsR_d;
  double *gF2Dr, *gF2Dr_d;
  double *gF2DrR, *gF2DrR_d;
  double *gF2Ds, *gF2Ds_d;
  double *gF2DsR, *gF2DsR_d;
  double *gFInterp0, *gFInterp0_d;
  double *gFInterp0R, *gFInterp0R_d;
  double *gFInterp1, *gFInterp1_d;
  double *gFInterp1R, *gFInterp1R_d;
  double *gFInterp2, *gFInterp2_d;
  double *gFInterp2R, *gFInterp2R_d;
  double *gInterp, *gInterp_d;

  double *invMass, *invMass_d;
  double *lift, *lift_d;
  double *mass, *mass_d;
  double *r, *r_d;
  double *s, *s_d;
  double *ones, *ones_d;
  #ifdef INS_CUDA
  cublasHandle_t handle;
  #endif
};

#endif
