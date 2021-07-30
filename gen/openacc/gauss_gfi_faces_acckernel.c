//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void gauss_gfi_faces_openacc( const int *edgeNum, const bool *rev,
                            double **gf0, double **gf1, double **gf2) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 3 + n;
        indR = m * 3 + n;
      } else {
        indL = m * 3 + n;
        indR = (3 - 1 - m) * 3 + n;
      }

      if(edgeL == 0) {
        if(edgeR == 0) {
          gf0[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp0_g[indL];
        } else if(edgeR == 1) {
          gf0[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp0_g[indL];
        } else {
          gf0[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp0_g[indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          gf1[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp1_g[indL];
        } else if(edgeR == 1) {
          gf1[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp1_g[indL];
        } else {
          gf1[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp1_g[indL];
        }
      } else {
        if(edgeR == 0) {
          gf2[0][indL] += gFInterp0_g[indR];
          gf0[1][indR] += gFInterp2_g[indL];
        } else if(edgeR == 1) {
          gf2[0][indL] += gFInterp1_g[indR];
          gf1[1][indR] += gFInterp2_g[indL];
        } else {
          gf2[0][indL] += gFInterp2_g[indR];
          gf2[1][indR] += gFInterp2_g[indL];
        }
      }
    }
  }
}

// host stub function
void op_par_loop_gauss_gfi_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6){

  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 9, "double", OP_INC);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 9, "double", OP_INC);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 9, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(11);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[11].name      = name;
  OP_kernels[11].count    += 1;

  int  ninds   = 3;
  int  inds[8] = {-1,-1,0,0,1,1,2,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_gfi_faces\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_11
    int part_size = OP_PART_SIZE_11;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    bool* data1 = (bool*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data4 = (double *)arg4.data_d;
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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4,data6)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        double* arg2_vec[] = {
           &data2[9 * map2idx],
           &data2[9 * map3idx]};
        double* arg4_vec[] = {
           &data4[9 * map2idx],
           &data4[9 * map3idx]};
        double* arg6_vec[] = {
           &data6[9 * map2idx],
           &data6[9 * map3idx]};

        gauss_gfi_faces_openacc(
          &data0[2 * n],
          &data1[1 * n],
          arg2_vec,
          arg4_vec,
          arg6_vec);
      }

    }
    OP_kernels[11].transfer  += Plan->transfer;
    OP_kernels[11].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[11].time     += wall_t2 - wall_t1;
}
