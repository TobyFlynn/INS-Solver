//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void viscosity_bc_openacc( const int *bedge_type, const int *bedgeNum,
                         const double *t, const double *x, const double *y,
                         double *exQ0, double *exQ1) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  if(*bedge_type == 0) {

    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 7; i++) {
      double y1 = y[exInd + i];
      exQ0[exInd + i] += pow(0.41, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (0.41 - y1);



    }
  } else if(*bedge_type == 1) {






  } else {






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
  op_arg arg6){

  double*arg2h = (double *)arg2.data;
  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(54);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[54].name      = name;
  OP_kernels[54].count    += 1;

  int  ninds   = 4;
  int  inds[7] = {-1,-1,-1,0,1,2,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: viscosity_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_54
    int part_size = OP_PART_SIZE_54;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg2_l = arg2h[0];

  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map3 = arg3.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map3,data0,data1,data3,data4,data5,data6)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map3idx;
        map3idx = map3[n + set_size1 * 0];


        viscosity_bc_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &arg2_l,
          &data3[21 * map3idx],
          &data4[21 * map3idx],
          &data5[21 * map3idx],
          &data6[21 * map3idx]);
      }

    }
    OP_kernels[54].transfer  += Plan->transfer;
    OP_kernels[54].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[54].time     += wall_t2 - wall_t1;
}
