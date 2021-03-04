//
// auto-generated by op2.py
//

//user function
#include "../kernels/poisson_test_error.h"

// host stub function
void op_par_loop_poisson_test_error(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg4h = (double *)arg4.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_error");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  // allocate and initialise arrays for global reduction
  double arg4_l[nthreads*64];
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg4_l[d+thr*64]=ZERO_double;
    }
  }

  if (set_size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        poisson_test_error(
          &((double*)arg0.data)[15*n],
          &((double*)arg1.data)[15*n],
          &((double*)arg2.data)[15*n],
          &((double*)arg3.data)[15*n],
          &arg4_l[64*omp_get_thread_num()]);
      }
    }
  }

  // combine reduction data
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg4h[d] += arg4_l[d+thr*64];
    }
  }
  op_mpi_reduce(&arg4,arg4h);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].time     += wall_t2 - wall_t1;
  OP_kernels[25].transfer += (float)set->size * arg0.size;
  OP_kernels[25].transfer += (float)set->size * arg1.size;
  OP_kernels[25].transfer += (float)set->size * arg2.size;
  OP_kernels[25].transfer += (float)set->size * arg3.size * 2.0f;
}
