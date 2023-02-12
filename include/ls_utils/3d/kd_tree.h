#ifndef __INS_KD_TREE_3D_H
#define __INS_KD_TREE_3D_H

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_mesh/dg_mesh_3d.h"
#include "poly_approx.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  double x_rot;
  double y_rot;
  double z_rot;
  double x;
  double y;
  double z;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  int axis;
};

class KDTree3D {
public:
  KDTree3D(const double *x, const double *y, const double *z, const int num,
           DGMesh3D *m, op_dat s);

  virtual void build_tree();
  KDCoord closest_point(double x, double y, double z);
  std::vector<PolyApprox3D> get_polys();

  std::vector<KDCoord> points;
  bool empty;
protected:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  double bb_sqr_dist(const int node_ind, const double x, const double y, const double z);
  void nearest_neighbour(double x, double y, double z, int current_ind, std::vector<KDCoord>::iterator &closest_pt, double &closest_distance);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<KDNode> nodes;

  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox3D> polys;
  static const int leaf_size = 10;
  DGMesh3D *mesh;
};

#endif