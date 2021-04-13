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

  virtual bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0) = 0;

  double getAverageConvergeIter();

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);

  // OP2 Dats
  op_dat bc_dat;
protected:
  void create_vec(Vec *v, int size = 15);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat, int size = 15);
  void store_vec(Vec *v, op_dat v_dat);
  void create_mat(Mat *m, int row, int col, int prealloc);
  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  int *dirichlet;
  int *neumann;

  bool massMat;
  double massFactor;

  int numberIter = 0;
  int solveCount = 0;
};

class Poisson_M : public Poisson {
public:
  Poisson_M(INSData *data, CubatureData *cubData, GaussData *gaussData) : Poisson(data, cubData, gaussData) {}
  ~Poisson_M();

  void createMatrix();
  void createMassMatrix();
  void createBCMatrix();

  bool solve(op_dat b_dat, op_dat x_dat, bool addMass = false, double factor = 0.0);

private:
  Mat pMat, pBCMat, pMMat;

  bool pMatInit = false;
  bool pMMatInit = false;
  bool pBCMatInit = false;
};

#endif
