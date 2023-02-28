#include "dg_compiler_defs.h"
#include "cblas.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

// global constants
extern int DG_CONSTANTS[6];
extern int FMASK[120];
extern double r_ynolds;
extern double mu0;
extern double mu1;
extern double rho0;
extern double rho1;

// header
#include "op_lib_cpp.h"

inline void _pmf_3d_mult_faces(const int *mapping, const int **order, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP *in, const DG_FP *in_x,
                              const DG_FP *in_y, const DG_FP *in_z, DG_FP *l_x,
                              DG_FP *l_y, DG_FP *l_z, DG_FP *out) {
  const int p = order[0][0];
  const int dg_npf = DG_CONSTANTS[(p - 1) * DG_NUM_CONSTANTS + 1];

  const DG_FP gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < dg_npf; j++) {
    const int fmaskL = mapping[j * 4 + 0];
    const int fmaskR = mapping[j * 4 + 1];
    const int indL = mapping[j * 4 + 2];
    const int indR = mapping[j * 4 + 3];
    const DG_FP diffL_u = in[fmaskL] - in[fmaskR];
    const DG_FP diffL_u_x = nx[0] * (in_x[fmaskR] + in_x[fmaskL]);
    const DG_FP diffL_u_y = ny[0] * (in_y[fmaskR] + in_y[fmaskL]);
    const DG_FP diffL_u_z = nz[0] * (in_z[fmaskR] + in_z[fmaskL]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    out[indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[indL] += nx[0] * l_tmpL;
    l_y[indL] += ny[0] * l_tmpL;
    l_z[indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[fmaskR] - in[fmaskL];
    const DG_FP diffR_u_x = nx[1] * (in_x[fmaskR] + in_x[fmaskL]);
    const DG_FP diffR_u_y = ny[1] * (in_y[fmaskR] + in_y[fmaskL]);
    const DG_FP diffR_u_z = nz[1] * (in_z[fmaskR] + in_z[fmaskL]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    out[indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[indR] += nx[1] * l_tmpR;
    l_y[indR] += ny[1] * l_tmpR;
    l_z[indR] += nz[1] * l_tmpR;
  }
}

// host stub function
void custom_cpu_kernel_pmf_3d_mult_faces(char const *name, op_set set,
  op_arg mapping,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16,
  op_arg arg18,
  op_arg arg20,
  op_arg arg22,
  op_arg arg24){

  int nargs = 26;
  op_arg args[26];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 20, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 20, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 20, "double", OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 20, "double", OP_READ);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 40, "double", OP_INC);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 40, "double", OP_INC);
  }

  arg22.idx = 0;
  args[22] = arg22;
  for ( int v=1; v<2; v++ ){
    args[22 + v] = op_arg_dat(arg22.dat, v, arg22.map, 40, "double", OP_INC);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<2; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, 40, "double", OP_INC);
  }

  int  ninds   = 9;
  int  inds[26] = {0,0,-1,-1,-1,-1,-1,-1,-1,-1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pmf_3d_mult_faces\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_38
    int part_size = OP_PART_SIZE_38;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
          int map0idx;
          int map1idx;
          map0idx = arg0.map_data[n * arg0.map->dim + 0];
          map1idx = arg0.map_data[n * arg0.map->dim + 1];

          const int* arg0_vec[] = {
             &((int*)arg0.data)[1 * map0idx],
             &((int*)arg0.data)[1 * map1idx]};

          _pmf_3d_mult_faces(
            &((int*)mapping.data)[4 * DG_NPF * n],
            arg0_vec,
            &((double*)arg5.data)[2 * n],
            &((double*)arg6.data)[2 * n],
            &((double*)arg7.data)[2 * n],
            &((double*)arg8.data)[2 * n],
            &((double*)arg9.data)[2 * n],
            (double *)arg10.data,
            (double *)arg12.data,
            (double *)arg14.data,
            (double *)arg16.data,
            (double *)arg18.data,
            (double *)arg20.data,
            (double *)arg22.data,
            (double *)arg24.data);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[38].transfer  += Plan->transfer;
    OP_kernels[38].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);
}
