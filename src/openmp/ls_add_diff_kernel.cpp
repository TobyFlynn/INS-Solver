//
// auto-generated by op2.py
//

//user function
#include "../kernels/ls_add_diff.h"

// host stub function
void op_par_loop_ls_add_diff(char const *name, op_set set,
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
  op_timing_realloc(68);
  OP_kernels[68].name      = name;
  OP_kernels[68].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_add_diff");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  if (set_size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        ls_add_diff(
          &((double*)arg0.data)[15*n],
          &((double*)arg1.data)[15*n],
          &((double*)arg2.data)[21*n],
          &((double*)arg3.data)[21*n],
          &((double*)arg4.data)[21*n],
          &((double*)arg5.data)[21*n]);
      }
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[68].time     += wall_t2 - wall_t1;
  OP_kernels[68].transfer += (float)set->size * arg0.size;
  OP_kernels[68].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[68].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[68].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[68].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[68].transfer += (float)set->size * arg5.size * 2.0f;
}
