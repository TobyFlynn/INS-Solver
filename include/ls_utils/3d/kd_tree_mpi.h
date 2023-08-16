#ifndef __INS_KD_TREE_MPI_3D_H
#define __INS_KD_TREE_MPI_3D_H

#include "dg_compiler_defs.h"
#include "kd_tree.h"
#include <vector>
#include "mpi.h"

class KDTree3DMPI : public KDTree3D {
public:
  KDTree3DMPI(DGMesh3D *m, const DG_FP alpha);

  void build_tree(const DG_FP *x, const DG_FP *y, const DG_FP *z, const int num,
                  op_dat s) override;
private:
  DG_FP min_dist_bb(const DG_FP *min_0, const DG_FP *max_0, const DG_FP *min_1,
                    const DG_FP *max_1);
  std::vector<int> ranks;

  MPI_Datatype packed_type;

  int rank, comm_size;
};

#endif
