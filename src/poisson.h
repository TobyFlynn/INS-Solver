#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"

extern double gFInterp0[7 * 15];
extern double gFInterp1[7 * 15];
extern double gFInterp2[7 * 15];
extern double gaussW[7];

class Poisson {
public:
  Poisson(INSData *nsData, CubatureData *cubData, GaussData *gaussData);
  ~Poisson();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);

  double getAverageConvergeIter();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  void createMatrix();
  void createMassMatrix();
  void createBCMatrix();

  // OP2 Dats
  op_dat bc_dat;
private:
  void create_vec(Vec *v, int size = 15);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat, int size = 15);
  void store_vec(Vec *v, op_dat v_dat);
  void create_mat(Mat *m, int row, int col, int prealloc);
  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  Mat pMat, pBCMat, pMMat;

  int *dirichlet;
  int *neumann;

  bool massMat;
  double massFactor;

  int numberIter = 0;
  int solveCount = 0;

  bool pMatInit = false;
  bool pMMatInit = false;
  bool pBCMatInit = false;
};

#endif
