//
// auto-generated by op2.py
//

void pressure_bc_omp4_kernel(
  int *data0,
  int dat0size,
  int *data1,
  int dat1size,
  double *arg2,
  int *arg3,
  int *map4,
  int map4size,
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  double arg2_l = *arg2;
  int arg3_l = *arg3;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: reynolds_ompkernel, FMASK_ompkernel[:9])\
    map(to:col_reord[0:set_size1],map4[0:map4size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size],data12[0:dat12size],data13[0:dat13size],data14[0:dat14size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map4idx;
    map4idx = map4[n_op + set_size1 * 0];

    //variable mapping
    const int *bedge_type = &data0[1*n_op];
    const int *bedgeNum = &data1[1*n_op];
    const double *t = &arg2_l;
    const int *problem = &arg3_l;
    const double *x = &data4[6 * map4idx];
    const double *y = &data5[6 * map4idx];
    const double *nx = &data6[9 * map4idx];
    const double *ny = &data7[9 * map4idx];
    const double *nu = &data8[6 * map4idx];
    const double *rho = &data9[6 * map4idx];
    const double *N0 = &data10[6 * map4idx];
    const double *N1 = &data11[6 * map4idx];
    const double *gradCurlVel0 = &data12[6 * map4idx];
    const double *gradCurlVel1 = &data13[6 * map4idx];
    double *dPdN = &data14[9 * map4idx];

    //inline function
    
    int exInd = *bedgeNum * 3;
    int *fmask = &FMASK_ompkernel[*bedgeNum * 3];

    const double PI = 3.141592653589793238463;

    if(*problem == 0) {
      if(*bedge_type == 0 || *bedge_type == 2 || *bedge_type == 3) {

        for(int i = 0; i < 3; i++) {
          int fInd = fmask[i];


          double res1 = -N0[fInd] - gradCurlVel1[fInd] / (reynolds_ompkernel * rho[fInd]);
          double res2 = -N1[fInd] + gradCurlVel0[fInd] / (reynolds_ompkernel * rho[fInd]);
          dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
        }
      }

      if(*bedge_type == 0) {

        for(int i = 0; i < 3; i++) {
          double y1 = y[fmask[i]];
          double bcdUndt = -pow(1.0, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0) * 6.0 * y1 * (1.0 - y1);
          dPdN[exInd + i] -= bcdUndt;
        }
      }
    } else {
      if(*bedge_type == 0) {

        for(int i = 0; i < 3; i++) {
          int fInd = fmask[i];
          double res1 = -N0[fInd] - nu[fInd] * gradCurlVel1[fInd];
          double res2 = -N1[fInd] + nu[fInd] * gradCurlVel0[fInd];
          dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;

          double y1 = y[fmask[i]];
          double x1 = x[fmask[i]];
          double nx1 = nx[exInd + i];
          double ny1 = ny[exInd + i];
          double bcdUndt = -nu[fInd] * 4.0 * PI * PI * (-nx1 * sin(2.0 * PI * y1) + ny1 * sin(2.0 * PI * x1))
                            * exp(-nu[fInd] * 4.0 * PI * PI * *t);
          dPdN[exInd + i] -= bcdUndt;
        }
      }
    }
    //end inline func
  }

  *arg2 = arg2_l;
  *arg3 = arg3_l;
}
