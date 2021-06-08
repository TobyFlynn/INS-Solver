//
// auto-generated by op2.py
//

void diff_bflux_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data6,
  int dat6size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size]) \
    map(to: gaussW_g_ompkernel[:7])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data6[0:dat6size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    map1idx = map1[n_op + set_size1 * 0];

    //variable mapping
    const int *bedgeNum = &data0[1*n_op];
    const double *sJ = &data1[21 * map1idx];
    const double *nx = &data2[21 * map1idx];
    const double *ny = &data3[21 * map1idx];
    const double *sigX = &data4[21 * map1idx];
    const double *sigY = &data4[21 * map1idx];
    double *flux = &data6[21 * map1idx];

    //inline function
    
    int exInd = 0;
    if(*bedgeNum == 1) exInd = 7;
    else if(*bedgeNum == 2) exInd = 2 * 7;

    for(int i = 0; i < 7; i++) {
      flux[exInd + i] += 0.5 * gaussW_g_ompkernel[i] * sJ[exInd + i] * (nx[exInd + i] * sigX[exInd + i] + ny[exInd + i] * sigY[exInd + i]);

    }
    //end inline func
  }

}
