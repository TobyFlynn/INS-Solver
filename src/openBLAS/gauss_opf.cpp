#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_gauss_opf(const int numCells, double *op, const double *gf,
                               const double *pD, const double *term0,
                               const double *term1, const double *term2) {
  for(int c = 0; c < numCells; c++) {
    double *op_c = op + c * 15 * 15;
    const double *gf_c = gf + c * 7 * 15;
    const double *pD_c = pD + c * 7 * 15;
    const double *term0_c = term0 + c * 7 * 15;
    const double *term1_c = term1 + c * 7 * 15;
    const double *term2_c = term2 + c * 7 * 15;

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 15, 15, 7, 1.0, term0_c, 15, gf_c, 15, 0.0, op_c, 15);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 15, 15, 7, 1.0, term1_c, 15, pD_c, 15, 1.0, op_c, 15);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 15, 15, 7, -1.0, term2_c, 15, gf_c, 15, 1.0, op_c, 15);
  }
}

void gauss_opf_blas(INSData *nsData, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg gauss_args[] = {
    // Face 0
    op_arg_dat(gaussData->OPf[0], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->pGF[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pD[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDx[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 1
    op_arg_dat(gaussData->OPf[1], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->pGF[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pD[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->mDy[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    // Face 2
    op_arg_dat(gaussData->OPf[2], -1, OP_ID, 15 * 15, "double", OP_WRITE),
    op_arg_dat(gaussData->pGF[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pD[2], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[0], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[1], -1, OP_ID, 7 * 15, "double", OP_READ),
    op_arg_dat(gaussData->pDx[2], -1, OP_ID, 7 * 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(nsData->cells, 18, gauss_args);

  openblas_gauss_opf(nsData->numCells, (double *)gaussData->OPf[0]->data,
                   (double *)gaussData->pGF[0]->data, (double *)gaussData->pD[0]->data, (double *)gaussData->mDx[0]->data,
                   (double *)gaussData->mDx[1]->data, (double *)gaussData->mDx[2]->data);

  openblas_gauss_opf(nsData->numCells, (double *)gaussData->OPf[1]->data,
                   (double *)gaussData->pGF[1]->data, (double *)gaussData->pD[1]->data, (double *)gaussData->mDy[0]->data,
                   (double *)gaussData->mDy[1]->data, (double *)gaussData->mDy[2]->data);

  openblas_gauss_opf(nsData->numCells, (double *)gaussData->OPf[2]->data,
                   (double *)gaussData->pGF[2]->data, (double *)gaussData->pD[2]->data, (double *)gaussData->pDx[0]->data,
                   (double *)gaussData->pDx[1]->data, (double *)gaussData->pDx[2]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(18, gauss_args);
}
