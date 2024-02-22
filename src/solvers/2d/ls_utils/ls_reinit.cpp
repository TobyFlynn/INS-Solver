#include "solvers/2d/ls_solver.h"

#include <vector>
#include <memory>
#include <random>

#ifdef INS_MPI
#include "ls_utils/2d/kd_tree_mpi.h"
#else
#include "ls_utils/2d/kd_tree.h"
#endif
#include "dg_global_constants/dg_global_constants_2d.h"
#include "timing.h"
#include "op2_utils.h"
#include "dg_constants/dg_constants_2d.h"

extern Timing *timer;
extern DGConstants *constants;

void rs2xy(DG_FP &sampleX, DG_FP &sampleY, const DG_FP *nodeX, const DG_FP *nodeY) {
  DG_FP r_ = sampleX;
  DG_FP s_ = sampleY;

  sampleX = 0.5 * (nodeX[1] * (1.0 + r_) + nodeX[2] * (1.0 + s_) - nodeX[0] * (s_ + r_));
  sampleY = 0.5 * (nodeY[1] * (1.0 + r_) + nodeY[2] * (1.0 + s_) - nodeY[0] * (s_ + r_));
}

void set_initial_sample_pts(const DG_FP *surface, const DG_FP *nodeX, const DG_FP *nodeY, DG_FP *sampleX, DG_FP *sampleY) {
  const int fmask_node0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];
  const DG_FP node0s = surface[fmask_node0];
  const DG_FP node1s = surface[fmask_node1];
  const DG_FP node2s = surface[fmask_node2];

  DG_FP end0x = NAN;
  DG_FP end0y = NAN;
  DG_FP end1x = NAN;
  DG_FP end1y = NAN;

  if((node0s > 0.0) != (node1s > 0.0)) {
    end0x = nodeX[0] - (node0s / (node0s - node1s)) * (nodeX[0] - nodeX[1]);
    end0y = nodeY[0] - (node0s / (node0s - node1s)) * (nodeY[0] - nodeY[1]);
  }

  if((node1s > 0.0) != (node2s > 0.0)) {
    if(isnan(end0x)) {
      end0x = nodeX[1] - (node1s / (node1s - node2s)) * (nodeX[1] - nodeX[2]);
      end0y = nodeY[1] - (node1s / (node1s - node2s)) * (nodeY[1] - nodeY[2]);
    } else {
      end1x = nodeX[1] - (node1s / (node1s - node2s)) * (nodeX[1] - nodeX[2]);
      end1y = nodeY[1] - (node1s / (node1s - node2s)) * (nodeY[1] - nodeY[2]);
    }
  }

  if((node2s > 0.0) != (node0s > 0.0)) {
    if(isnan(end0x)) {
      end0x = nodeX[2] - (node2s / (node2s - node0s)) * (nodeX[2] - nodeX[0]);
      end0y = nodeY[2] - (node2s / (node2s - node0s)) * (nodeY[2] - nodeY[0]);
    } else {
      end1x = nodeX[2] - (node2s / (node2s - node0s)) * (nodeX[2] - nodeX[0]);
      end1y = nodeY[2] - (node2s / (node2s - node0s)) * (nodeY[2] - nodeY[0]);
    }
  }

  if(isnan(end0x) || isnan(end0y) || isnan(end1x) || isnan(end1y)) {
    const DG_FP *r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
    const DG_FP *s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;

    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sampleX[i] = r_ptr[i] * 0.75;
      sampleY[i] = s_ptr[i] * 0.75;
      rs2xy(sampleX[i], sampleY[i], nodeX, nodeY);
    }
    return;
  }

  DG_FP dist = sqrt((end1x - end0x) * (end1x - end0x) + (end1y - end0y) * (end1y - end0y));
  DG_FP dist_per_sample = dist / (LS_SAMPLE_NP - 1.0);

  DG_FP incrementx = ((end1x - end0x) / dist) * dist_per_sample;
  DG_FP incrementy = ((end1y - end0y) / dist) * dist_per_sample;

  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sampleX[i] = end0x + incrementx * i;
    sampleY[i] = end0y + incrementy * i;
  }
}

DG_FP in_tri_sign(const DG_FP x, const DG_FP y, const DG_FP v0x, const DG_FP v0y, const DG_FP v1x, const DG_FP v1y) {
  return (x - v1x) * (v0y - v1y) - (v0x - v1x) * (y - v1y);
}

bool in_tri(const DG_FP ptX, const DG_FP ptY, const DG_FP *nodeX, const DG_FP *nodeY) {
  bool has_neg, has_pos;

  DG_FP d1 = in_tri_sign(ptX, ptY, nodeX[0], nodeY[0], nodeX[1], nodeY[1]);
  DG_FP d2 = in_tri_sign(ptX, ptY, nodeX[1], nodeY[1], nodeX[2], nodeY[2]);
  DG_FP d3 = in_tri_sign(ptX, ptY, nodeX[2], nodeY[2], nodeX[0], nodeY[0]);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);
}

bool simplified_newton(DG_FP &pt_x, DG_FP &pt_y, PolyApprox &pol, const DG_FP tol) {
  bool converged = false;
  for(int step = 0; step < 100; step++) {
    DG_FP surf = pol.val_at(pt_x, pt_y);
    DG_FP dsdx, dsdy;
    pol.grad_at(pt_x, pt_y, dsdx, dsdy);

    DG_FP sqrnorm = dsdx * dsdx + dsdy * dsdy;
    if(sqrnorm > 1e-14) {
      dsdx *= surf / sqrnorm;
      dsdy *= surf / sqrnorm;
    }

    pt_x -= dsdx;
    pt_y -= dsdy;

    // Check convergence
    if(dsdx * dsdx + dsdy * dsdy < tol) {
      converged = true;
      break;
    }
  }
  return converged;
}

void LevelSetSolver2D::sampleInterface(op_dat sampleX, op_dat sampleY, std::vector<PolyApprox> &polys, 
                                       std::map<int,int> &cell2polyMap, std::set<int> &cellInds) {
  // Setup random number generator for later
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  const DG_FP tol = fmax(1e-18, 1e-4 * h * h);
  // const DG_FP tol = 1e-16;

  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  const DG_FP *nodeX_ptr = getOP2PtrHost(mesh->nodeX, OP_READ);
  const DG_FP *nodeY_ptr = getOP2PtrHost(mesh->nodeY, OP_READ);
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  DG_FP *sampleX_ptr = getOP2PtrHost(sampleX, OP_WRITE);
  DG_FP *sampleY_ptr = getOP2PtrHost(sampleY, OP_WRITE);

  // Reset all sample points
  #pragma omp parallel for
  for(int cell = 0; cell < mesh->cells->size; cell++) {
    DG_FP *sampleX_c = sampleX_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleY_c = sampleY_ptr + cell * LS_SAMPLE_NP;
    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sampleX_c[i] = NAN;
      sampleY_c[i] = NAN;
    }
  }

  // Set sample points
  std::vector<int> cell_inds_vec(cellInds.begin(), cellInds.end());
  #pragma omp parallel for
  for(int cell = 0; cell < cell_inds_vec.size(); cell++) {
    const int cell_ind = cell_inds_vec[cell];
    const DG_FP *surface = s_ptr + cell_ind * DG_NP;
    const DG_FP *nodeX = nodeX_ptr + cell_ind * 3;
    const DG_FP *nodeY = nodeY_ptr + cell_ind * 3;
    DG_FP *sampleX = sampleX_ptr + cell_ind * LS_SAMPLE_NP;
    DG_FP *sampleY = sampleY_ptr + cell_ind * LS_SAMPLE_NP;

    // Set initial position of sample points
    set_initial_sample_pts(surface, nodeX, nodeY, sampleX, sampleY);

    // Simplified Newton method
    int poly_ind = cell2polyMap.at(cell_ind);
    PolyApprox pol = polys[poly_ind];
    DG_FP off_x, off_y;
    pol.get_offsets(off_x, off_y);

    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sampleX[i] -= off_x;
      sampleY[i] -= off_y;
    }

    const DG_FP _nodeX[] = {nodeX[0] - off_x, nodeX[1] - off_x, nodeX[2] - off_x};
    const DG_FP _nodeY[] = {nodeY[0] - off_y, nodeY[1] - off_y, nodeY[2] - off_y};

    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      bool converged = false;
      DG_FP start_x = sampleX[i];
      DG_FP start_y = sampleY[i];
      if(!isnan(sampleX[i]))
        converged = simplified_newton(sampleX[i], sampleY[i], pol, tol);
      DG_FP dist_travelled = (start_x - sampleX[i]) * (start_x - sampleX[i]) + (start_y - sampleY[i]) * (start_y - sampleY[i]);

      // || !in_tri(sampleX[i], sampleY[i], _nodeX, _nodeY)
      // || dist_travelled > 1.5 * 1.5 * h * h
      if(!converged || dist_travelled > 1.5 * 1.5 * h * h) {
        // Try randomly placing start and rerunning 10 times
        bool rerun_converged = false;
        bool rerun_in_bounds = false;
        int rerun_counter = 0;
        while((!rerun_converged || !rerun_in_bounds) && rerun_counter < 10) {
          // Random start point
          sampleX[i] = dis(gen);
          sampleY[i] = dis(gen);
          rs2xy(sampleX[i], sampleY[i], _nodeX, _nodeY);
          while(!in_tri(sampleX[i], sampleY[i], _nodeX, _nodeY)) {
            sampleX[i] = dis(gen);
            sampleY[i] = dis(gen);
            rs2xy(sampleX[i], sampleY[i], _nodeX, _nodeY);
          }
          // Rerun
          DG_FP start_x = sampleX[i];
          DG_FP start_y = sampleY[i];
          rerun_converged = simplified_newton(sampleX[i], sampleY[i], pol, tol);
          DG_FP dist_travelled = (start_x - sampleX[i]) * (start_x - sampleX[i]) + (start_y - sampleY[i]) * (start_y - sampleY[i]);
          rerun_in_bounds = dist_travelled > 1.5 * 1.5 * h * h;
          rerun_counter++;
        }

        if(!rerun_converged || !rerun_in_bounds) {
          sampleX[i] = NAN;
          sampleY[i] = NAN;
        }
      }
    }

    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sampleX[i] += off_x;
      sampleY[i] += off_y;
    }
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(s, OP_READ, s_ptr);
  releaseOP2PtrHost(mesh->nodeX, OP_READ, nodeX_ptr);
  releaseOP2PtrHost(mesh->nodeY, OP_READ, nodeY_ptr);
  releaseOP2PtrHost(sampleX, OP_WRITE, sampleX_ptr);
  releaseOP2PtrHost(sampleY, OP_WRITE, sampleY_ptr);
}