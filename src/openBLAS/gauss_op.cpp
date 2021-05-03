#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_gauss_op(const int numCells, double *op, const double *mD,
                              const double *term0, const double *term1,
                              const double *term2, int face) {
  double *gFInterp;
  if(face == 0) {
    gFInterp = constants->gFInterp0;
  } else if (face == 1) {
    gFInterp = constants->gFInterp1;
  } else {
    gFInterp = constants->gFInterp2;
  }

  for(int c = 0; c < numCells; c++) {
    double *op_c = op + c * 15 * 15;
    const double *mD_c = mD + c * 7 * 15;
    const double *term0_c = term0 + c * 7 * 15;
    const double *term1_c = term1 + c * 7 * 15;
    const double *term2_c = term2 + c * 7 * 15;

    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 15, 15, 7, 1.0, term0_c, 7, gFInterp, 15, 0.0, op_c, 15);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 15, 15, 7, -1.0, term1_c, 7, mD_c, 15, 1.0, op_c, 15);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 15, 15, 7, -1.0, term2_c, 7, gFInterp, 15, 1.0, op_c, 15);
  }
}

void gauss_op_blas(INSData *nsData, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg gauss_args[] = {
    // Face 0
    op_arg_dat(gaussData->OP[0], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 1
    op_arg_dat(gaussData->OP[1], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 2
    op_arg_dat(gaussData->OP[2], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->mD[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[2], -1, OP_ID, 7 * 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(nsData->cells, 15, gauss_args);

  openblas_gauss_op(nsData->numCells, (double *)gaussData->OP[0]->data,
                   (double *)gaussData->mD[0]->data, (double *)gaussData->mDx[0]->data,
                   (double *)gaussData->mDx[1]->data, (double *)gaussData->mDx[2]->data, 0);

  openblas_gauss_op(nsData->numCells, (double *)gaussData->OP[1]->data,
                   (double *)gaussData->mD[1]->data, (double *)gaussData->mDy[0]->data,
                   (double *)gaussData->mDy[1]->data, (double *)gaussData->mDy[2]->data, 1);

  openblas_gauss_op(nsData->numCells, (double *)gaussData->OP[2]->data,
                   (double *)gaussData->mD[2]->data, (double *)gaussData->pDx[0]->data,
                   (double *)gaussData->pDx[1]->data, (double *)gaussData->pDx[2]->data, 2);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(15, gauss_args);
}
