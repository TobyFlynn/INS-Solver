#ifndef __INS_KD_TREE_MPI_3D_H
#define __INS_KD_TREE_MPI_3D_H

#include "kd_tree.h"

class KDTree3DMPI : public KDTree3D {
public:
  KDTree3DMPI(const double *x, const double *y, const double *z, const int num,
              DGMesh3D *m, op_dat s);

  void build_tree() override;
private:
};

#endif