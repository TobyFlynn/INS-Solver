#include "solvers/3d/ls_solver.h"

#include "op_seq.h"

#include "ls_utils/3d/kd_tree.h"
#ifdef INS_MPI
#include "ls_utils/3d/kd_tree_mpi.h"
#endif
#include "utils.h"
#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "dg_dat_pool.h"

#include <iostream>
#include <fstream>

#include "timing.h"

extern Timing *timer;
extern DGDatPool *dg_dat_pool;

LevelSetSolver3D::LevelSetSolver3D(DGMesh3D *m) {
  mesh = m;
  resuming = false;
  advectionSolver = new AdvectionSolver3D(m);

  s = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_solver_s");

  h = 0;
  reinit_count = 0;
}

LevelSetSolver3D::LevelSetSolver3D(DGMesh3D *m, const std::string &filename) {
  mesh = m;
  resuming = true;
  advectionSolver = new AdvectionSolver3D(m);

  s = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ls_solver_s");

  h = 0;
  reinit_count = 0;
}

LevelSetSolver3D::~LevelSetSolver3D() {
  delete advectionSolver;
}

void LevelSetSolver3D::setBCTypes(op_dat bc) {
  advectionSolver->set_bc_types(bc);
}

void LevelSetSolver3D::init() {
  if(!resuming) {
    op_par_loop(init_surface_3d, "init_surface_3d", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  // alpha = 2.0 * h / DG_ORDER;
  // order_width = 2.0 * h;
  // epsilon = h / DG_ORDER;
  alpha = 12.0 * h;
  // order_width = 12.0 * h;
  // epsilon = h;
  // reinit_width = 10.0 * h;
  // reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  // numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  reinitLS();
}

void LevelSetSolver3D::getRhoMu(op_dat rho, op_dat mu) {
  timer->startTimer("LevelSetSolver3D - getRhoMu");
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver3D - getRhoMu");
}

void LevelSetSolver3D::getNormalsCurvature(op_dat nx, op_dat ny, op_dat nz, op_dat curv) {
  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  timer->startTimer("LevelSetSolver3D - getNormalsCurvature");
  mesh->grad(s, nx, ny, nz);
  mesh->div(nx, ny, nz, curv);
  timer->endTimer("LevelSetSolver3D - getNormalsCurvature");
}

void LevelSetSolver3D::step(op_dat u, op_dat v, op_dat w, const DG_FP dt, const int num_steps) {
  timer->startTimer("LevelSetSolver3D - step");
  advectionSolver->set_dt(dt);
  for(int i = 0; i < num_steps; i++)
    advectionSolver->step(s, u, v, w);

  op_par_loop(ls_post_advec, "ls_post_advec", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  if(reinit_count > 49) {
    reinitLS();
    reinit_count = 0;
    op_par_loop(ls_post_reinit, "ls_post_reinit", mesh->cells,
                op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }
  reinit_count++;
  timer->endTimer("LevelSetSolver3D - step");
}

bool newtoncp_gepp(arma::mat &A, arma::vec &b) {
  for (int i = 0; i < 4; ++i) {
    int j = i;
    for (int k = i + 1; k < 4; ++k)
      if (std::abs(A(k,i)) > std::abs(A(j,i)))
        j = k;
    if (j != i) {
      for (int k = 0; k < 4; ++k)
        std::swap(A(i,k), A(j,k));
      std::swap(b(i), b(j));
    }

    if (std::abs(A(i,i)) < 1.0e4*std::numeric_limits<DG_FP>::epsilon())
      return false;

    DG_FP fac = 1.0 / A(i,i);
    for (int j = i + 1; j < 4; ++j)
      A(j,i) *= fac;

    for (int j = i + 1; j < 4; ++j) {
      for (int k = i + 1; k < 4; ++k)
        A(j,k) -= A(j,i)*A(i,k);
      b(j) -= A(j,i)*b(i);
    }
  }

  for (int i = 4 - 1; i >= 0; --i) {
    DG_FP sum = 0.0;
    for (int j = i + 1; j < 4; ++j)
      sum += A(i,j)*b(j);
    b(i) = (b(i) - sum) / A(i,i);
  }

  return true;
}

bool newton_kernel(DG_FP &closest_pt_x, DG_FP &closest_pt_y, DG_FP &closest_pt_z,
                   const DG_FP node_x, const DG_FP node_y, const DG_FP node_z,
                   PolyApprox3D &p, const DG_FP h) {
  DG_FP lambda = 0.0;
  bool converged = false;
  DG_FP pt_x = closest_pt_x;
  DG_FP pt_y = closest_pt_y;
  DG_FP pt_z = closest_pt_z;
  DG_FP init_x = closest_pt_x;
  DG_FP init_y = closest_pt_y;
  DG_FP init_z = closest_pt_z;

  for(int step = 0; step < 100; step++) {
    DG_FP pt_x_old = pt_x;
    DG_FP pt_y_old = pt_y;
    DG_FP pt_z_old = pt_z;
    // Evaluate surface and gradient at current guess
    DG_FP surface = p.val_at(pt_x, pt_y, pt_z);
    DG_FP surface_dx, surface_dy, surface_dz;
    p.grad_at(pt_x, pt_y, pt_z, surface_dx, surface_dy, surface_dz);
    // Evaluate Hessian
    DG_FP hessian[6];
    p.hessian_at(pt_x, pt_y, pt_z, hessian[0], hessian[1], hessian[2],
                 hessian[3], hessian[4], hessian[5]);

    // Check if |nabla(surface)| = 0, if so then return
    DG_FP gradsqrnorm = surface_dx * surface_dx + surface_dy * surface_dy + surface_dz * surface_dz;
    if(gradsqrnorm < 1e-14)
      break;

    // Init lambda at first step
    if(step == 0)
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy + (node_z - pt_z) * surface_dz) / gradsqrnorm;

    // Gradient of functional
    arma::vec gradf(4);
    gradf(0) = pt_x - node_x + lambda * surface_dx;
    gradf(1) = pt_y - node_y + lambda * surface_dy;
    gradf(2) = pt_z - node_z + lambda * surface_dz;
    gradf(3) = surface;

    // Calculate Hessian of functional
    arma::mat hessianf(4, 4);
    hessianf(0, 0) = 1.0 + lambda * hessian[0];
    hessianf(0, 1) = lambda * hessian[3]; hessianf(1, 0) = hessianf(0, 1);
    hessianf(0, 2) = lambda * hessian[4]; hessianf(2, 0) = hessianf(0, 2);
    hessianf(0, 3) = surface_dx; hessianf(3, 0) = hessianf(0, 3);

    hessianf(1, 1) = 1.0 + lambda * hessian[1];
    hessianf(1, 2) = lambda * hessian[5]; hessianf(2, 1) = hessianf(1, 2);
    hessianf(1, 3) = surface_dy; hessianf(3, 1) = hessianf(1, 3);

    hessianf(2, 2) = 1.0 + lambda * hessian[2];
    hessianf(2, 3) = surface_dz; hessianf(3, 2) = hessianf(2, 3);

    hessianf(3, 3) = 0.0;

    if(!newtoncp_gepp(hessianf, gradf)) {
      DG_FP delta1_x = (surface / gradsqrnorm) * surface_dx;
      DG_FP delta1_y = (surface / gradsqrnorm) * surface_dy;
      DG_FP delta1_z = (surface / gradsqrnorm) * surface_dz;
      lambda = ((node_x - pt_x) * surface_dx + (node_y - pt_y) * surface_dy + (node_z - pt_z) * surface_dz) / gradsqrnorm;
      DG_FP delta2_x = pt_x - node_x + lambda * surface_dx;
      DG_FP delta2_y = pt_y - node_y + lambda * surface_dy;
      DG_FP delta2_z = pt_z - node_z + lambda * surface_dz;
      DG_FP msqr = delta2_x * delta2_x + delta2_y * delta2_y + delta2_z * delta2_z;
      if(msqr > 0.1 * h * 0.1 * h) {
        delta2_x *= 0.1 * h / sqrt(msqr);
        delta2_y *= 0.1 * h / sqrt(msqr);
        delta2_z *= 0.1 * h / sqrt(msqr);
      }
      pt_x -= delta1_x + delta2_x;
      pt_y -= delta1_y + delta2_y;
      pt_z -= delta1_z + delta2_z;
    } else {
      arma::vec ans = gradf;

      // Clamp update
      DG_FP msqr = ans(0) * ans(0) + ans(1) * ans(1) + ans(2) * ans(2);
      if(msqr > h * 0.5 * h * 0.5)
        ans = ans * 0.5 * h / sqrt(msqr);

      // Update guess
      pt_x -= ans(0);
      pt_y -= ans(1);
      pt_z -= ans(2);
      lambda -= ans(3);
    }

    // Gone outside the element, return
    // if((init_x - pt_x) * (init_x - pt_x) + (init_y - pt_y) * (init_y - pt_y) > h * h) {
    //   pt_x = pt_x_old;
    //   pt_y = pt_y_old;
    //   pt_z = pt_z_old;
    //   break;
    // }

    // Converged, no more steps required
    if(sqrt((pt_x_old - pt_x) * (pt_x_old - pt_x) + (pt_y_old - pt_y) * (pt_y_old - pt_y) + (pt_z_old - pt_z) * (pt_z_old - pt_z)) < 1e-12) {
      converged = true;
      break;
    }
  }

  closest_pt_x = pt_x;
  closest_pt_y = pt_y;
  closest_pt_z = pt_z;

  return converged;
}

void newton_method(const int numPts, DG_FP *closest_x, DG_FP *closest_y, DG_FP *closest_z,
                   const DG_FP *x, const DG_FP *y, const DG_FP *z, int *poly_ind,
                   std::vector<PolyApprox3D> &polys, DG_FP *s, const DG_FP h) {
  int numNonConv = 0;
  int numReinit = 0;

  #pragma omp parallel for
  for(int i = 0; i < numPts; i++) {
    int start_ind = (i / DG_NP) * DG_NP;
    bool reinit = false;
    for(int j = 0; j < DG_NP; j++) {
      if(fabs(s[start_ind + j]) < 0.02) {
        reinit = true;
      }
    }

    if(reinit) {
      bool tmp = newton_kernel(closest_x[i], closest_y[i], closest_z[i], x[i], y[i], z[i], polys[poly_ind[i]], h);
      if(tmp) {
        bool negative = s[i] < 0.0;
        s[i] = (closest_x[i] - x[i]) * (closest_x[i] - x[i]) + (closest_y[i] - y[i]) * (closest_y[i] - y[i]) + (closest_z[i] - z[i]) * (closest_z[i] - z[i]);
        s[i] = sqrt(s[i]);
        if(negative) s[i] *= -1.0;
      }
      if(!tmp) {
        #pragma omp atomic
        numNonConv++;
      }
      #pragma omp atomic
      numReinit++;
    }
  }

  if(numNonConv != 0)
    std::cout << numNonConv << " non-converged points out of " << numReinit << " points reinitialised" << std::endl;
}

void LevelSetSolver3D::reinitLS() {
  #ifdef DG_OP2_SOA
  throw std::runtime_error("reinitLS not implemented for SoA");
  #endif
  timer->startTimer("LevelSetSolver3D - reinitLS");
  DGTempDat tmp_sampleX = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleY = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  DGTempDat tmp_sampleZ = dg_dat_pool->requestTempDatCells(LS_SAMPLE_NP);
  timer->startTimer("LevelSetSolver3D - sample interface");
  sampleInterface(tmp_sampleX.dat, tmp_sampleY.dat, tmp_sampleZ.dat);
  timer->endTimer("LevelSetSolver3D - sample interface");

  timer->startTimer("LevelSetSolver3D - build KD-Tree");
  const DG_FP *sample_pts_x = getOP2PtrHost(tmp_sampleX.dat, OP_READ);
  const DG_FP *sample_pts_y = getOP2PtrHost(tmp_sampleY.dat, OP_READ);
  const DG_FP *sample_pts_z = getOP2PtrHost(tmp_sampleZ.dat, OP_READ);

  #ifdef INS_MPI
  KDTree3DMPI kdtree(sample_pts_x, sample_pts_y, sample_pts_z, LS_SAMPLE_NP * mesh->cells->size, mesh, s);
  kdtree.build_tree();
  #else
  KDTree3D kdtree(sample_pts_x, sample_pts_y, sample_pts_z, LS_SAMPLE_NP * mesh->cells->size, mesh, s);
  kdtree.build_tree();
  #endif

  releaseOP2PtrHost(tmp_sampleX.dat, OP_READ, sample_pts_x);
  releaseOP2PtrHost(tmp_sampleY.dat, OP_READ, sample_pts_y);
  releaseOP2PtrHost(tmp_sampleZ.dat, OP_READ, sample_pts_z);

  dg_dat_pool->releaseTempDatCells(tmp_sampleX);
  dg_dat_pool->releaseTempDatCells(tmp_sampleY);
  dg_dat_pool->releaseTempDatCells(tmp_sampleZ);
  timer->endTimer("LevelSetSolver3D - build KD-Tree");

  timer->startTimer("LevelSetSolver3D - query KD-Tree");
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHost(mesh->z, OP_READ);

  DG_FP *closest_x = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_y = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *closest_z = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  int *poly_ind     = (int *)calloc(DG_NP * mesh->cells->size, sizeof(int));
  std::vector<PolyApprox3D> polys;

  if(!kdtree.empty) {
    #pragma omp parallel for
    for(int i = 0; i < DG_NP * mesh->cells->size; i++) {
      // Get closest sample point
      KDCoord tmp = kdtree.closest_point(x_ptr[i], y_ptr[i], z_ptr[i]);
      closest_x[i] = tmp.x;
      closest_y[i] = tmp.y;
      closest_z[i] = tmp.z;
      poly_ind[i]  = tmp.poly;
    }

    polys = kdtree.get_polys();
  }
  timer->endTimer("LevelSetSolver3D - query KD-Tree");

  if(h == 0.0) {
    op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
                op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
                op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
    h = 1.0 / h;
  }

  timer->startTimer("LevelSetSolver3D - newton method");
  DG_FP *surface_ptr = getOP2PtrHost(s, OP_RW);
  // Newton method
  if(!kdtree.empty) {
    newton_method(DG_NP * mesh->cells->size, closest_x, closest_y, closest_z,
                  x_ptr, y_ptr, z_ptr, poly_ind, polys, surface_ptr, h);
  }
  releaseOP2PtrHost(s, OP_RW, surface_ptr);
  timer->endTimer("LevelSetSolver3D - newton method");

  free(closest_x);
  free(closest_y);
  free(closest_z);
  free(poly_ind);

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);
  timer->endTimer("LevelSetSolver3D - reinitLS");
}
