//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf2_opf_openacc( const double *tol, const int *edgeNum, const double *gop0L,
                            const double *gop1L, const double *gop2L, const double *gopf0L,
                            const double *gopf1L, const double *gopf2L, double *op2L,
                            double *op1L,
                            const double *gop0R, const double *gop1R, const double *gop2R,
                            const double *gopf0R, const double *gopf1R, const double *gopf2R,
                            double *op2R, double *op1R) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(edgeL == 0) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop0L[colInd];
        if(fabs(val) > *tol) {
          op1L[ind] += val;
        }
        val = -0.5 * gopf0L[colInd];
        if(fabs(val) > *tol) {
          op2L[ind] += val;
        }
      }
    }
  } else if(edgeL == 1) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop1L[colInd];
        if(fabs(val) > *tol) {
          op1L[ind] += val;
        }
        val = -0.5 * gopf1L[colInd];
        if(fabs(val) > *tol) {
          op2L[ind] += val;
        }
      }
    }
  } else {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop2L[colInd];
        if(fabs(val) > *tol) {
          op1L[ind] += val;
        }
        val = -0.5 * gopf2L[colInd];
        if(fabs(val) > *tol) {
          op2L[ind] += val;
        }
      }
    }
  }

  if(edgeR == 0) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop0R[colInd];
        if(fabs(val) > *tol) {
          op1R[ind] += val;
        }
        val = -0.5 * gopf0R[colInd];
        if(fabs(val) > *tol) {
          op2R[ind] += val;
        }
      }
    }
  } else if(edgeR == 1) {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop1R[colInd];
        if(fabs(val) > *tol) {
          op1R[ind] += val;
        }
        val = -0.5 * gopf1R[colInd];
        if(fabs(val) > *tol) {
          op2R[ind] += val;
        }
      }
    }
  } else {
    for(int m = 0; m < 15; m++) {
      for(int n = 0; n < 15; n++) {
        int ind = m * 15 + n;
        int colInd = n * 15 + m;
        double val = 0.5 * gop2R[colInd];
        if(fabs(val) > *tol) {
          op1R[ind] += val;
        }
        val = -0.5 * gopf2R[colInd];
        if(fabs(val) > *tol) {
          op2R[ind] += val;
        }
      }
    }
  }
}

// host stub function
void op_par_loop_poisson_mf2_opf(char const *name, op_set set,
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
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17){

  double*arg0h = (double *)arg0.data;
  int nargs = 18;
  op_arg args[18];

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
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(26);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[26].name      = name;
  OP_kernels[26].count    += 1;

  int  ninds   = 7;
  int  inds[18] = {-1,-1,0,1,2,3,4,5,-1,6,0,1,2,3,4,5,-1,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf2_opf\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_26
    int part_size = OP_PART_SIZE_26;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data1 = (int*)arg1.data_d;
    double* data8 = (double*)arg8.data_d;
    double* data16 = (double*)arg16.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data9 = (double *)arg9.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data1,data8,data16,data2,data3,data4,data5,data6,data7,data9)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map10idx;
        map2idx = map2[n + set_size1 * 0];
        map10idx = map2[n + set_size1 * 1];


        poisson_mf2_opf_openacc(
          &arg0_l,
          &data1[2 * n],
          &data2[225 * map2idx],
          &data3[225 * map2idx],
          &data4[225 * map2idx],
          &data5[225 * map2idx],
          &data6[225 * map2idx],
          &data7[225 * map2idx],
          &data8[225 * n],
          &data9[225 * map2idx],
          &data2[225 * map10idx],
          &data3[225 * map10idx],
          &data4[225 * map10idx],
          &data5[225 * map10idx],
          &data6[225 * map10idx],
          &data7[225 * map10idx],
          &data16[225 * n],
          &data9[225 * map10idx]);
      }

    }
    OP_kernels[26].transfer  += Plan->transfer;
    OP_kernels[26].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[26].time     += wall_t2 - wall_t1;
}
