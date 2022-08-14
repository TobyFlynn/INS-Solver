#include "op_seq.h"

#include <memory>
#include <iostream>
#include <stdexcept>

#include "dg_utils.h"

double *getOP2PtrDevice(op_dat dat, op_access acc) {
  throw std::runtime_error("\ngetOP2PtrDevice not implemented for CPU\n");
}

void releaseOP2PtrDevice(op_dat dat, op_access acc, const double *ptr) {
  throw std::runtime_error("\releaseOP2PtrDevice not implemented for CPU\n");
}

double *getOP2PtrHost(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges(dat->set, 1, args);
  op_mpi_wait_all(1, args);
  return (double *)dat->data;
}

void releaseOP2PtrHost(op_dat dat, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };

  op_mpi_set_dirtybit(1, args);

  ptr = nullptr;
}

double *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc) {
  throw std::runtime_error("\ngetOP2PtrDevice not implemented for CPU\n");
}

void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const double *ptr) {
  throw std::runtime_error("\releaseOP2PtrDevice not implemented for CPU\n");
}

double *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges(map->from, 2, args);
  op_mpi_wait_all(2, args);
  return (double *)dat->data;
}

void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };

  op_mpi_set_dirtybit(2, args);

  ptr = nullptr;
}

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_SUB_CELLS, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * DG_SUB_CELLS * sizeof(double));
  memcpy(res, dat->data, dat->set->size * DG_SUB_CELLS * sizeof(double));
  op_mpi_set_dirtybit(1, args);
  return res;
}

bool is_point_in_cell(const double x, const double y, const double *cellX, const double *cellY) {
  double ABx = cellX[1] - cellX[0];
  double ABy = cellY[1] - cellY[0];

  double APx = x - cellX[0];
  double APy = y - cellY[0];

  double BCx = cellX[2] - cellX[1];
  double BCy = cellY[2] - cellY[1];

  double BPx = x - cellX[1];
  double BPy = y - cellY[1];

  double CAx = cellX[0] - cellX[2];
  double CAy = cellY[0] - cellY[2];

  double CPx = x - cellX[2];
  double CPy = y - cellY[2];

  double AB_AP = ABx * APy - ABy * APx;
  double BC_BP = BCx * BPy - BCy * BPx;
  double CA_CP = CAx * CPy - CAy * CPx;

  bool zero0 = AB_AP == 0.0 || AB_AP == -0.0;
  bool zero1 = BC_BP == 0.0 || BC_BP == -0.0;
  bool zero2 = CA_CP == 0.0 || CA_CP == -0.0;

  if(zero0 || zero1 || zero2) {
    // One zero means on an edge of the triangle
    // Two zeros means on a vertex of the triangle
    if(zero0 && !zero1 && !zero2) {
      return (BC_BP > 0.0 && CA_CP > 0.0) || (BC_BP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && zero1 && !zero2) {
      return (AB_AP > 0.0 && CA_CP > 0.0) || (AB_AP < 0.0 && CA_CP < 0.0);
    } else if(!zero0 && !zero1 && zero2) {
      return (AB_AP > 0.0 && BC_BP > 0.0) || (AB_AP < 0.0 && BC_BP < 0.0);
    } else if(zero0 && zero1 && !zero2) {
      return true;
    } else if(zero0 && !zero1 && zero2) {
      return true;
    } else if(!zero0 && zero1 && zero2) {
      return true;
    } else {
      return false;
    }
  }

  if(AB_AP > 0.0 && BC_BP > 0.0 && CA_CP > 0.0) {
    return true;
  } else if(AB_AP < 0.0 && BC_BP < 0.0 && CA_CP < 0.0) {
    return true;
  } else {
    return false;
  }
}
