//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_edges_openacc( const double *uL, const double *opL, double *rhsL,
                          const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < 10; m++) {
    int ind = m * 10;
    for(int n = 0; n < 10; n++) {
      rhsL[m] += opL[ind + n] * uR[n];
      rhsR[m] += opR[ind + n] * uL[n];
    }
  }
}

// host stub function
void op_par_loop_poisson_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(44);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[44].name      = name;
  OP_kernels[44].count    += 1;

  int  ninds   = 2;
  int  inds[6] = {0,-1,1,0,-1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_edges\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_44
    int part_size = OP_PART_SIZE_44;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map0 = arg0.map_data_d;

    double* data1 = (double*)arg1.data_d;
    double* data4 = (double*)arg4.data_d;
    double *data0 = (double *)arg0.data_d;
    double *data2 = (double *)arg2.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map0,data1,data4,data0,data2)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map0idx;
        int map3idx;
        map0idx = map0[n + set_size1 * 0];
        map3idx = map0[n + set_size1 * 1];


        poisson_edges_openacc(
          &data0[10 * map0idx],
          &data1[100 * n],
          &data2[10 * map0idx],
          &data0[10 * map3idx],
          &data4[100 * n],
          &data2[10 * map3idx]);
      }

    }
    OP_kernels[44].transfer  += Plan->transfer;
    OP_kernels[44].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[44].time     += wall_t2 - wall_t1;
}
