//
// auto-generated by op2.py
//

void poisson_set_rhs_bc_omp4_kernel(
  int *data0,
  int dat0size,
  int *data1,
  int dat1size,
  int *arg2,
  int *arg3,
  int *map4,
  int map4size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  int arg2_l = *arg2;
  int arg3_l = *arg3;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: FMASK_ompkernel[:15])\
    map(to:col_reord[0:set_size1],map4[0:map4size],data4[0:dat4size],data5[0:dat5size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map4idx;
    map4idx = map4[n_op + set_size1 * 0];

    //variable mapping
    const int *bedge_type = &data0[1*n_op];
    const int *bedgeNum = &data1[1*n_op];
    const int *dirichlet0 = &arg2_l;
    const int *dirichlet1 = &arg3_l;
    const double *tau = &data4[15 * map4idx];
    double *bcTau = &data5[15 * map4idx];

    //inline function
    
    int exInd = 0;
    if(*bedgeNum == 1) {
      exInd = 5;
    } else if(*bedgeNum == 2) {
      exInd = 2 * 5;
    }

    int *fmask;

    if(*bedgeNum == 0) {
      fmask = FMASK_ompkernel;
    } else if(*bedgeNum == 1) {
      fmask = &FMASK_ompkernel[5];
    } else {
      fmask = &FMASK_ompkernel[2 * 5];
    }

    if(*bedge_type == *dirichlet0 || *bedge_type == *dirichlet1) {
      for(int i = 0; i < 5; i++) {
        bcTau[exInd + i] += tau[exInd + i];
      }
    }
    //end inline func
  }

  *arg2 = arg2_l;
  *arg3 = arg3_l;
}
