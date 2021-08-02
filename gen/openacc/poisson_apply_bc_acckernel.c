//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_apply_bc_openacc( const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * 6;

  for(int m = 0; m < 10; m++) {
    int ind = m * 6;
    for(int n = 0; n < 6; n++) {
      rhs[m] += op[ind + n] * bc[exInd + n];
    }
  }
}

// host stub function
void op_par_loop_poisson_apply_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;

  int  ninds   = 2;
  int  inds[4] = {-1,-1,0,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_apply_bc\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_15
    int part_size = OP_PART_SIZE_15;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data3)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        map2idx = map2[n + set_size1 * 0];


        poisson_apply_bc_openacc(
          &data0[1 * n],
          &data1[60 * n],
          &data2[18 * map2idx],
          &data3[10 * map2idx]);
      }

    }
    OP_kernels[15].transfer  += Plan->transfer;
    OP_kernels[15].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
}
