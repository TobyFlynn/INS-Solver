#include "constants.h"

Constants::Constants() {
  // Cubature constants
  cubDr  = cubDr_g;
  cubDs  = cubDs_g;
  cubV   = cubV_g;
  cubVDr = cubVDr_g;
  cubVDs = cubVDs_g;
  cubW   = cubW_g;
  // Grad constants
  Dr  = Dr_g;
  Drw = Drw_g;
  Ds  = Ds_g;
  Dsw = Dsw_g;
  // Gauss constants
  gaussW     = gaussW_g;
  gF0Dr      = gF0Dr_g;
  gF0DrR     = gF0DrR_g;
  gF0Ds      = gF0Ds_g;
  gF0DsR     = gF0DsR_g;
  gF1Dr      = gF1Dr_g;
  gF1DrR     = gF1DrR_g;
  gF1Ds      = gF1Ds_g;
  gF1DsR     = gF1DsR_g;
  gF2Dr      = gF2Dr_g;
  gF2DrR     = gF2DrR_g;
  gF2Ds      = gF2Ds_g;
  gF2DsR     = gF2DsR_g;
  gFInterp0  = gFInterp0_g;
  gFInterp0R = gFInterp0R_g;
  gFInterp1  = gFInterp1_g;
  gFInterp1R = gFInterp1R_g;
  gFInterp2  = gFInterp2_g;
  gFInterp2R = gFInterp2R_g;
  gInterp    = gInterp_g;
  // Other constants
  invMass = invMass_g;
  lift    = lift_g;
  mass    = mass_g;
  r       = r_g;
  s       = s_g;
  ones    = ones_g;
}

Constants::~Constants() {

}

double* Constants::get_ptr(Constant_Matrix mat) {
  switch(mat) {
    case CUB_DR:
      return cubDr;
    case CUB_DS:
      return cubDs;
    case CUB_V:
      return cubV;
    case CUB_VDR:
      return cubVDr;
    case CUB_VDS:
      return cubVDs;
    case CUB_W:
      return cubW;
    case DR:
      return Dr;
    case DRW:
      return Drw;
    case DS:
      return Ds;
    case DSW:
      return Dsw;
    case GAUSS_W:
      return gaussW;
    case GAUSS_F0DR:
      return gF0Dr;
    case GAUSS_F0DR_R:
      return gF0DrR;
    case GAUSS_F0DS:
      return gF0Ds;
    case GAUSS_F0DS_R:
      return gF0DsR;
    case GAUSS_F1DR:
      return gF1Dr;
    case GAUSS_F1DR_R:
      return gF1DrR;
    case GAUSS_F1DS:
      return gF1Ds;
    case GAUSS_F1DS_R:
      return gF1DsR;
    case GAUSS_F2DR:
      return gF2Dr;
    case GAUSS_F2DR_R:
      return gF2DrR;
    case GAUSS_F2DS:
      return gF2Ds;
    case GAUSS_F2DS_R:
      return gF2DsR;
    case GAUSS_FINTERP0:
      return gFInterp0;
    case GAUSS_FINTERP0_R:
      return gFInterp0R;
    case GAUSS_FINTERP1:
      return gFInterp1;
    case GAUSS_FINTERP1_R:
      return gFInterp1R;
    case GAUSS_FINTERP2:
      return gFInterp2;
    case GAUSS_FINTERP2_R:
      return gFInterp2R;
    case GAUSS_INTERP:
      return gInterp;
    case INV_MASS:
      return invMass;
    case LIFT:
      return lift;
    case MASS:
      return mass;
    case R:
      return r;
    case S:
      return s;
    case ONES:
      return ones;
  }
  return nullptr;
}
