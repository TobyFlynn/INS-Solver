//
// auto-generated by op2.py
//

void poisson_mf_edges_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  int *map2,
  int map2size,
  double *data2,
  int dat2size,
  double *data4,
  int dat4size,
  double *data6,
  int dat6size,
  double *data8,
  int dat8size,
  double *data10,
  int dat10size,
  double *data12,
  int dat12size,
  double *data14,
  int dat14size,
  double *data16,
  int dat16size,
  double *data18,
  int dat18size,
  double *data20,
  int dat20size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: gaussW_g_ompkernel[:7])\
    map(to:col_reord[0:set_size1],map2[0:map2size],data2[0:dat2size],data4[0:dat4size],data6[0:dat6size],data8[0:dat8size],data10[0:dat10size],data12[0:dat12size],data14[0:dat14size],data16[0:dat16size],data18[0:dat18size],data20[0:dat20size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map2idx;
    int map3idx;
    map2idx = map2[n_op + set_size1 * 0];
    map3idx = map2[n_op + set_size1 * 1];

    const double* arg2_vec[] = {
       &data2[21 * map2idx],
       &data2[21 * map3idx]};
    const double* arg4_vec[] = {
       &data4[21 * map2idx],
       &data4[21 * map3idx]};
    const double* arg6_vec[] = {
       &data6[21 * map2idx],
       &data6[21 * map3idx]};
    const double* arg8_vec[] = {
       &data8[3 * map2idx],
       &data8[3 * map3idx]};
    const double* arg10_vec[] = {
       &data10[21 * map2idx],
       &data10[21 * map3idx]};
    const double* arg12_vec[] = {
       &data12[21 * map2idx],
       &data12[21 * map3idx]};
    const double* arg14_vec[] = {
       &data14[21 * map2idx],
       &data14[21 * map3idx]};
    double* arg16_vec[] = {
       &data16[21 * map2idx],
       &data16[21 * map3idx]};
    double* arg18_vec[] = {
       &data18[21 * map2idx],
       &data18[21 * map3idx]};
    double* arg20_vec[] = {
       &data20[21 * map2idx],
       &data20[21 * map3idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    const double **sJ = arg2_vec;
    const double **nx = arg4_vec;
    const double **ny = arg6_vec;
    const double **tau = arg8_vec;
    const double **u = arg10_vec;
    const double **dudx = arg12_vec;
    const double **dudy = arg14_vec;
    double **fluxX = arg16_vec;
    double **fluxY = arg18_vec;
    double **flux = arg20_vec;

    //inline function
    
    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;

    int exIndL = 0;
    if(edgeL == 1) exIndL = 7;
    else if(edgeL == 2) exIndL = 2 * 7;

    int exIndR = 0;
    if(edgeR == 1) exIndR = 7;
    else if(edgeR == 2) exIndR = 2 * 7;

    for(int i = 0; i < 7; i++) {
      int rInd;
      int lInd = exIndL + i;
      if(reverse) {
        rInd = exIndR + 7 - i - 1;
      } else {
        rInd = exIndR + i;
      }

      double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
      tmp *= gaussW_g_ompkernel[i] * sJ[0][lInd];
      fluxX[0][lInd] += nx[0][lInd] * tmp;
      fluxY[0][lInd] += ny[0][lInd] * tmp;
      tmp = nx[0][lInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
      tmp += ny[0][lInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
      tmp -= tau[0][edgeL] * (u[0][lInd] - u[1][rInd]) / 2.0;
      tmp *= gaussW_g_ompkernel[i] * sJ[0][lInd];
      flux[0][lInd] += tmp;
    }

    for(int i = 0; i < 7; i++) {
      int lInd;
      int rInd = exIndR + i;
      if(reverse) {
        lInd = exIndL + 7 - i - 1;
      } else {
        lInd = exIndL + i;
      }

      double tmp = (u[0][lInd] + u[1][rInd]) / 2.0;
      tmp *= gaussW_g_ompkernel[i] * sJ[1][rInd];
      fluxX[1][rInd] += nx[1][rInd] * tmp;
      fluxY[1][rInd] += ny[1][rInd] * tmp;
      tmp = nx[1][rInd] * ((dudx[0][lInd] + dudx[1][rInd]) / 2.0);
      tmp += ny[1][rInd] * ((dudy[0][lInd] + dudy[1][rInd]) / 2.0);
      tmp -= tau[1][edgeR] * (u[1][rInd] - u[0][lInd]) / 2.0;
      tmp *= gaussW_g_ompkernel[i] * sJ[1][rInd];
      flux[1][rInd] += tmp;
    }
    //end inline func
  }

}
