//
// auto-generated by op2.py
//

//user function
#include "../kernels/calc_h.h"

// host stub function
void op_par_loop_calc_h(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  double*arg2h = (double *)arg2.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(60);
  OP_kernels[60].name      = name;
  OP_kernels[60].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  calc_h");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  // allocate and initialise arrays for global reduction
  double arg2_l[nthreads*64];
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg2_l[d+thr*64]=arg2h[d];
    }
  }

  if (set_size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        calc_h(
          &((double*)arg0.data)[3*n],
          &((double*)arg1.data)[3*n],
          &arg2_l[64*omp_get_thread_num()]);
      }
    }
  }

  // combine reduction data
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg2h[d]  = MIN(arg2h[d],arg2_l[d+thr*64]);
    }
  }
  op_mpi_reduce(&arg2,arg2h);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[60].time     += wall_t2 - wall_t1;
  OP_kernels[60].transfer += (float)set->size * arg0.size;
  OP_kernels[60].transfer += (float)set->size * arg1.size;
}
