inline void poisson_mf2_faces(const int *edgeNum, const double *uL,
                              const double *op0L, const double *op1L,
                              const double *op2L, double *rhsL,
                              const double *uR, const double *op0R, const double *op1R,
                              const double *op2R, double *rhsR) {
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
/*
  const double *opL;
  if(edgeL == 0) {
    opL = op0[0];
  } else if(edgeL == 1) {
    opL = op1[0];
  } else {
    opL = op2[0];
  }

  const double *opR;
  if(edgeR == 0) {
    opR = op0[1];
  } else if(edgeR == 1) {
    opR = op1[1];
  } else {
    opR = op2[1];
  }
*/
  if(edgeL == 0) {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op0L[ind + n] * uR[n];
      }
      rhsL[m] += val;
    }
  } else if(edgeL == 1) {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op1L[ind + n] * uR[n];
      }
      rhsL[m] += val;
    }
  } else {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op2L[ind + n] * uR[n];
      }
      rhsL[m] += val;
    }
  }

  if(edgeR == 0) {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op0R[ind + n] * uL[n];
      }
      rhsR[m] += val;
    }
  } else if(edgeR == 1) {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op1R[ind + n] * uL[n];
      }
      rhsR[m] += val;
    }
  } else {
    for(int m = 0; m < 15; m++) {
      int ind = m * 15;
      double val = 0.0;
      for(int n = 0; n < 15; n++) {
        val += op2R[ind + n] * uL[n];
      }
      rhsR[m] += val;
    }
  }

/*
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opL[ind + n] * u[1][n];
    }
    rhs[0][m] += val;
  }

  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {
      val += opR[ind + n] * u[0][n];
    }
    rhs[1][m] += val;
  }
*/
}
