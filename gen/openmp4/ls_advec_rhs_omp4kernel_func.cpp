//
// auto-generated by op2.py
//

void ls_advec_rhs_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  double *data11,
  int dat11size,
  double *data12,
  int dat12size,
  double *data13,
  int dat13size,
  double *data14,
  int dat14size,
  double *data15,
  int dat15size,
  double *data16,
  int dat16size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size],data12[0:dat12size],data13[0:dat13size],data14[0:dat14size],data15[0:dat15size],data16[0:dat16size]) \
    map(to: FMASK_ompkernel[:12])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *dFdr = &data0[10*n_op];
    const double *dFds = &data1[10*n_op];
    const double *dGdr = &data2[10*n_op];
    const double *dGds = &data3[10*n_op];
    const double *rx = &data4[10*n_op];
    const double *ry = &data5[10*n_op];
    const double *sx = &data6[10*n_op];
    const double *sy = &data7[10*n_op];
    const double *q = &data8[10*n_op];
    const double *exQ = &data9[12*n_op];
    const double *u = &data10[10*n_op];
    const double *v = &data11[10*n_op];
    const double *fscale = &data12[12*n_op];
    const double *nx = &data13[12*n_op];
    const double *ny = &data14[12*n_op];
    double *nFlux = &data15[12*n_op];
    double *output = &data16[10*n_op];

    //inline function
    
    for(int i = 0; i < 10; i++) {
      output[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
    }

    double mQ[3 * 4];
    double mF[3 * 4];
    double mG[3 * 4];
    for(int i = 0; i < 3 * 4; i++) {
      int ind = FMASK_ompkernel[i];
      mQ[i] = q[ind];
      mF[i] = u[ind] * q[ind];
      mG[i] = v[ind] * q[ind];
    }

    double pF[3 * 4];
    double pG[3 * 4];
    for(int i = 0; i < 3 * 4; i++) {
      int ind = FMASK_ompkernel[i];
      pF[i]  = u[ind] * exQ[i];
      pG[i]  = v[ind] * exQ[i];
    }

    for(int i = 0; i < 3 * 4; i++) {
      int ind = FMASK_ompkernel[i];



      nFlux[i] = (nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] + exQ[i]);
      nFlux[i] += fabs(nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] - exQ[i]);
      nFlux[i] *= 0.5 * fscale[i];
    }
    //end inline func
  }

}
