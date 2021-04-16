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
  LIFT    = LIFT_g;
  MASS    = MASS_g;
  r       = r_g;
  s       = s_g;
  ones    = ones_g;
}

Constants::~Constants() {

}
