#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"
#include "petscvec.h"
#include "petscksp.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];
extern double Drw[15 * 15];
extern double Dsw[15 * 15];
extern double LIFT[15 * 15];
extern double MASS[15 * 15];
extern int FMASK[15];

class Poisson {
public:
  Poisson(INSData *nsData);
  ~Poisson();

  void rhs(const double *u, double *rhs);
  void solve(op_dat b_dat, op_dat x_dat, bool method, bool addMass = false, double factor = 0.0);

  void setDirichletBCs(int *d, op_dat d_dat);
  void setNeumannBCs(int *n);
  // OP2 Dats
  op_dat pTau, pExRHS[2], pU, pDu, pFluxXu, pFluxYu, pDuDx, pDuDy, pFluxQ, pDivQ, pRHS;
  op_dat dBC;
private:
  void copy_u(const double *u);
  void copy_rhs(double *rhs);
  void create_vec(Vec *v);
  void destroy_vec(Vec *v);
  void load_vec(Vec *v, op_dat v_dat);
  void store_vec(Vec *v, op_dat v_dat);
  INSData *data;
  // Pointers to private memory
  double *pTau_data;
  double *pExRHS_data[2];
  double *pU_data;
  double *pDu_data;
  double *pFluxXu_data;
  double *pFluxYu_data;
  double *pDuDx_data;
  double *pDuDy_data;
  double *pFluxQ_data;
  double *pDivQ_data;
  double *pRHS_data;

  int *dirichlet;
  int *neumann;

  bool massMat;
  double massFactor;
};

PetscErrorCode matAMult(Mat A, Vec x, Vec y);

#endif
