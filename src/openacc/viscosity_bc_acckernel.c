//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void viscosity_bc_openacc( const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, const double *nx, const double *ny,
                         const double *nu, double *exQ0, double *exQ1) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0) {

      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        exQ0[exInd + i] += pow(1.0, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (1.0 - y1);
      }
    } else if(*bedge_type == 1) {






    } else {

    }
  } else {
    if(*bedge_type == 0) {

      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        exQ0[exInd + i] += -sin(2.0 * PI * y1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += sin(2.0 * PI * x1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
      }
    }

    if(*bedge_type == 1) {

      for(int i = 0; i < 7; i++) {
        double y1  = y[exInd + i];
        double x1  = x[exInd + i];
        double ny1 = ny[exInd + i];
        double nx1 = nx[exInd + i];
        exQ0[exInd + i] += ny1 * 2.0 * PI * (-cos(2.0 * PI * y1)) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += nx1 * 2.0 * PI * cos(2.0 * PI * x1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
      }
    }
  }
}

// host stub function
void op_par_loop_viscosity_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10){

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(36);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[36].name      = name;
  OP_kernels[36].count    += 1;

  int  ninds   = 7;
  int  inds[11] = {-1,-1,-1,-1,0,1,2,3,4,5,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: viscosity_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_36
    int part_size = OP_PART_SIZE_36;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg2_l = arg2h[0];
  int arg3_l = arg3h[0];

  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map4 = arg4.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data9 = (double *)arg9.data_d;
    double *data10 = (double *)arg10.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map4,data0,data1,data4,data5,data6,data7,data8,data9,data10)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map4idx;
        map4idx = map4[n + set_size1 * 0];


        viscosity_bc_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &arg3_l,
          &data4[21 * map4idx],
          &data5[21 * map4idx],
          &data6[21 * map4idx],
          &data7[21 * map4idx],
          &data8[21 * map4idx],
          &data9[21 * map4idx],
          &data10[21 * map4idx]);
      }

    }
    OP_kernels[36].transfer  += Plan->transfer;
    OP_kernels[36].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[36].time     += wall_t2 - wall_t1;
}
