// Include VTK stuff
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellIterator.h>
#include <vtkIdList.h>
// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <getopt.h>

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"file", required_argument, 0, 0},
  {"bc", required_argument, 0, 0},
  {0,    0,                  0,  0}
};

// Structs for converting mesh
struct Point2D {
  double x;
  double y;
};

struct Cell {
  int points[3];
};

struct Edge {
  int points[2];
  int cells[2];
  int num[2];
};

int getBoundaryEdgeNum(const string &type, double x0, double y0, double x1, double y1) {
  if(type == "cylinder") {
    if(x0 == 0.0 && x1 == 0.0) {
      // Inflow
      return 0;
    } else if(x0 == x1 && x0 > 5.0) {
      // Outflow
      return 1;
    } else if(x0 > 0.1 && x1 > 0.1 && x0 < 1.0 && x1 < 1.0
              && y0 > 0.1 && y1 > 0.1 && y0 < 0.9 && y1 < 0.9) {
      // Cylinder Wall
      return 2;
    } else {
      // Top/Bottom Wall
      return 3;
    }
  } else if(type == "airfoil") {
    if(x0 == x1 && x0 < -4.9) {
      // Inflow
      return 0;
    } else if(x0 == x1 && x0 > 4.9) {
      // Outflow
      return 1;
    } else if(x0 > -1.0 && x1 > -1.0 && x0 < 2.0 && x1 < 2.0
              && y0 > -1.0 && y1 > -1.0 && y0 < 1.0 && y1 < 1.0) {
      // Airfoil Wall
      return 2;
    } else {
      cout << "TOP/BOTTOM Wall" << endl;
      // Top/Bottom Wall
      return 2;
    }
  } else if(type == "poisson-test-0") {
    if(y0 == y1 && y0 > 0.5) {
      // Neumann BC y = 1
      return 1;
    } else if(y0 == y1 && y0 < 0.5) {
      // Neumann BC y = 0
      return 1;
    } else if(x0 < 0.5){
      // Dirichlet BC x = 0
      return 1;
    } else {
      // Dirichlet BC x = 1
      return 0;
    }
  } else if(type == "poisson-test-1") {
    if(y0 == y1 && y0 > 0.5) {
      // Neumann BC y = 1
      return 2;
    } else if(y0 == y1 && y0 < 0.5) {
      // Neumann BC y = 0
      return 3;
    } else if(x0 < 0.5){
      // Dirichlet BC x = 0
      return 0;
    } else {
      // Dirichlet BC x = 1
      return 1;
    }
  } else if(type == "vortex") {
    if(y0 == -0.5 && y1 == -0.5 && x0 <= 1e-11 && x1 <= 1e-11) {
      // Outflow
      return 1;
    } else if(y0 == -0.5 && y1 == -0.5 && x0 >= -1e-11 && x1 >= -1e-11) {
      // Inflow
      return 0;
    } else if(y0 == 0.5 && y1 == 0.5 && x0 <= 1e-11 && x1 <= 1e-11) {
      // Inflow
      return 0;
    } else if(y0 == 0.5 && y1 == 0.5 && x0 >= -1e-11 && x1 >= -1e-11) {
      // Outflow
      return 1;
    } else if(x0 == -0.5 && x1 == -0.5 && y0 <= 1e-11 && y1 <= 1e-11) {
      // Inflow
      return 0;
    } else if(x0 == -0.5 && x1 == -0.5 && y0 >= -1e-11 && y1 >= -1e-11) {
      // Outflow
      return 1;
    } else if(x0 == 0.5 && x1 == 0.5 && y0 <= 1e-11 && y1 <= 1e-11) {
      // Outflow
      return 1;
    } else if(x0 == 0.5 && x1 == 0.5 && y0 >= -1e-11 && y1 >= -1e-11) {
      // Inflow
      return 0;
    } else {
      cerr << "***ERROR*** Boundary edge not categorised (vortex)" << endl;
      cerr << "   " << x0 << "," << y0 << " " << x1 << "," << y1 <<  endl;
      return -1;
    }
  } else {
    cerr << "***ERROR*** Unrecognised boundary type specified" << endl;
  }
  return -1;
}

void getBCs(const string &type, int *bc_data) {
  if(type == "cylinder") {
    // Pressure Dirichlet
    bc_data[0] = 1; bc_data[1] = -1; bc_data[2] = -1;
    // Pressure Neumann
    bc_data[3] = 0; bc_data[4] = 2; bc_data[5] = 3;
    // Viscosity Dirichlet
    bc_data[6] = 0; bc_data[7] = 2; bc_data[8] = 3;
    // Viscosity Neumann
    bc_data[9] = 1; bc_data[10] = -1; bc_data[11] = -1;
  } else if(type == "airfoil") {
    // // Pressure Dirichlet
    // bc_data[0] = 1; bc_data[1] = 3; bc_data[2] = -1;
    // // Pressure Neumann
    // bc_data[3] = 0; bc_data[4] = 2; bc_data[5] = -1;
    // // Viscosity Dirichlet
    // bc_data[6] = 0; bc_data[7] = 2; bc_data[8] = -1;
    // // Viscosity Neumann
    // bc_data[9] = 1; bc_data[10] = 3; bc_data[11] = -1;
    // Pressure Dirichlet
    bc_data[0] = 1; bc_data[1] = -1; bc_data[2] = -1;
    // Pressure Neumann
    bc_data[3] = 0; bc_data[4] = 2; bc_data[5] = 3;
    // Viscosity Dirichlet
    bc_data[6] = 0; bc_data[7] = 2; bc_data[8] = 3;
    // Viscosity Neumann
    bc_data[9] = 1; bc_data[10] = -1; bc_data[11] = -1;
  } else if(type == "poisson-test-0") {
    // N/A
  } else if(type == "poisson-test-1") {
    // N/A
  } else if(type == "vortex") {
    // Pressure Dirichlet
    bc_data[0] = 1; bc_data[1] = -1; bc_data[2] = -1;
    // Pressure Neumann
    bc_data[3] = 0; bc_data[4] = -1; bc_data[5] = -1;
    // Viscosity Dirichlet
    bc_data[6] = 0; bc_data[7] = -1; bc_data[8] = -1;
    // Viscosity Neumann
    bc_data[9] = 1; bc_data[10] = -1; bc_data[11] = -1;
  } else {
    cerr << "***ERROR*** Unrecognised boundary type specified" << endl;
  }
}

int main(int argc, char **argv) {
  string fileName = "naca0012.vtk";
  string bcType = "";
  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"file") == 0) fileName = optarg;
    if(strcmp((char*)options[opt_index].name,"bc") == 0) bcType = optarg;
  }

  // Read in VTK file generated by Gmsh
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName (fileName.c_str());
  reader->Update();
  grid = reader->GetOutput();

  map<int,unique_ptr<Cell>> cellMap;
  map<pair<int,int>,unique_ptr<Edge>> internalEdgeMap;

  vector<int> pointIds;
  vector<double> x;
  vector<double> y;

  // Iterate over cells
  vtkSmartPointer<vtkCellIterator> cellIterator = grid->NewCellIterator();
  while(!cellIterator->IsDoneWithTraversal()) {
    // Check that this is a cell (not a line or point)
    // only an issue due to how Gmsh generates the VTK file
    if(cellIterator->GetNumberOfPoints() == 3 && cellIterator->GetCellType() == VTK_TRIANGLE) {
      vtkSmartPointer<vtkIdList> ids = cellIterator->GetPointIds();
      int newIds[3];
      for(int i = 0; i < 3; i++) {
        auto it = find(pointIds.begin(), pointIds.end(), ids->GetId(i));
        if(it != pointIds.end()) {
          newIds[i] = distance(pointIds.begin(), it);
        } else {
          pointIds.push_back(ids->GetId(i));
          newIds[i] = pointIds.size() - 1;
          double coords[3];
          grid->GetPoint(ids->GetId(i), coords);
          x.push_back(coords[0]);
          y.push_back(coords[1]);
        }
      }

      // Add cell to map
      unique_ptr<Cell> cell = make_unique<Cell>();
      cell->points[0] = newIds[0] + 1;
      cell->points[1] = newIds[2] + 1;
      cell->points[2] = newIds[1] + 1;
      cellMap.insert(pair<int,unique_ptr<Cell>>(cellIterator->GetCellId(), move(cell)));
    }
    // Go to next cell
    cellIterator->GoToNextCell();
  }

  vector<cgsize_t> elements;
  for(auto const &elem : cellMap) {
    elements.push_back(elem.second->points[0]);
    elements.push_back(elem.second->points[1]);
    elements.push_back(elem.second->points[2]);
  }

  cout << "Number of points: " << x.size() << endl;
  cout << "VTK Number of points: " << grid->GetNumberOfPoints() << endl;
  cout << "Number of cell: " << elements.size() / 3 << endl;

  // Add edges to edge map if not already contained in mapping
  // If already added then update cell field
  // Try both combinations
  for(int i = 0; i < elements.size() / 3; i++) {
    int ind = i * 3;
    // Check that points are anticlockwise
    int p_1 = elements[ind];
    int p_2 = elements[ind + 1];
    int p_3 = elements[ind + 2];
    double val = (x[p_2 - 1] - x[p_1 - 1]) * (y[p_3 - 1] - y[p_2 - 1]) - (y[p_2 - 1] - y[p_1 - 1]) * (x[p_3 - 1] - x[p_2 - 1]);
    if(val < 0) {
      cout << "Switching points" << endl;
      cout << "  Old val: " << val << endl;
      elements[ind] = p_3;
      elements[ind + 2] = p_1;
      p_1 = elements[ind];
      p_2 = elements[ind + 1];
      p_3 = elements[ind + 2];
      val = (x[p_2 - 1] - x[p_1 - 1]) * (y[p_3 - 1] - y[p_2 - 1]) - (y[p_2 - 1] - y[p_1 - 1]) * (x[p_3 - 1] - x[p_2 - 1]);
      cout << "  New val: " << val << endl;
    }
    for(int j = 0; j < 3; j++) {
      int p1, p2;
      if(j == 0) {
        p1 = elements[ind];
        p2 = elements[ind + 1];
      } else if(j == 1) {
        p1 = elements[ind + 1];
        p2 = elements[ind + 2];
      } else {
        p1 = elements[ind + 2];
        p2 = elements[ind];
      }

      pair<int,int> key;
      if(p1 < p2) {
        key = make_pair(p1, p2);
      } else {
        key = make_pair(p2, p1);
      }

      if(internalEdgeMap.count(key) == 0) {
        unique_ptr<Edge> edge = make_unique<Edge>();
        edge->points[0] = key.first;
        edge->points[1] = key.second;
        edge->cells[0] = i;
        edge->cells[1] = -1;
        edge->num[0] = j;
        edge->num[1] = -1;
        internalEdgeMap.insert(make_pair(key, move(edge)));
      } else {
          if(internalEdgeMap.at(key)->cells[1] != -1) {
            cout << "ERROR in edge mapping: " << endl;
            cout << "  Old values: " << internalEdgeMap.at(key)->cells[0] << " " << internalEdgeMap.at(key)->cells[1] << endl;
            cout << "  New Value: " << i << endl;
            cout << "  Edges: " << internalEdgeMap.at(key)->points[0] << " " << internalEdgeMap.at(key)->points[1] << endl;
            cout << "  Key: " << key.first << " " << key.second << endl;
          }
          if(internalEdgeMap.at(key)->points[0] != key.first || internalEdgeMap.at(key)->points[1] != key.second) {
            cout << "ERROR in edge mapping: " << endl;
            cout << "  Prev Edge Nodes: " << internalEdgeMap.at(key)->points[0] << " " << internalEdgeMap.at(key)->points[1] << endl;
            cout << "  Current Nodes: " << p1 << " " << p2 << endl;
            cout << "  Key: " << key.first << " " << key.second << endl;
          }
          internalEdgeMap.at(key)->cells[1] = i;
          internalEdgeMap.at(key)->num[1] = j;
      }
    }
  }

  vector<int> edges;
  vector<int> boundaryEdges;
  for(auto const &edge : internalEdgeMap) {
    if(edge.second->cells[1] == -1) {
      boundaryEdges.push_back(edge.second->points[0]);
      boundaryEdges.push_back(edge.second->points[1]);
      boundaryEdges.push_back(edge.second->cells[0]);
      boundaryEdges.push_back(edge.second->num[0]);
      double x0 = x[edge.second->points[0] - 1];
      double y0 = y[edge.second->points[0] - 1];
      double x1 = x[edge.second->points[1] - 1];
      double y1 = y[edge.second->points[1] - 1];
      int bType = getBoundaryEdgeNum(bcType, x0, y0, x1, y1);
      boundaryEdges.push_back(bType);
    } else {
      if(edge.second->points[0] == edge.second->points[1])
        cout << "***** ERROR: Edge with identical points *****" << endl;
      if(edge.second->cells[0] == edge.second->cells[1])
        cout << "***** ERROR: Edge with identical cells *****" << endl;
      edges.push_back(edge.second->points[0]);
      edges.push_back(edge.second->points[1]);
      edges.push_back(edge.second->cells[0]);
      edges.push_back(edge.second->cells[1]);
      edges.push_back(edge.second->num[0]);
      edges.push_back(edge.second->num[1]);
    }
  }

  int bc_data[12];
  getBCs(bcType, bc_data);

  // Write out CGNS file
  string outFileName = "";
  for(int i = 0; i < fileName.length(); i++) {
    if(fileName.at(i) == '.')
      break;
    outFileName += fileName[i];
  }
  outFileName += ".cgns";
  int file;
  if (cg_open(outFileName.c_str(), CG_MODE_WRITE, &file)) {
    cg_error_exit();
  }

  // Create base
  string baseName = "Base";
  int cellDim = 2;
  int physicalDim = 2;
  int baseIndex;
  cg_base_write(file, baseName.c_str(), cellDim, physicalDim, &baseIndex);
  // Create zone
  string zoneName = "Zone1";
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = x.size();
  // Number of cells
  sizes[1] = cellMap.size();
  // Number of boundary vertices (zero if elements not sorted)
  sizes[2] = 0;
  int zoneIndex;
  cg_zone_write(file, baseIndex, zoneName.c_str(), sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndex;
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateX", x.data(), &coordIndex);
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateY", y.data(), &coordIndex);
  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];
  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                   start, end, 0, elements.data(), &sectionIndex);
  // Write edges
  // {p1, p2, c1, c2, num1, num2}
  int numEdges = edges.size() / 6;
  cout << "Number of edges: " << numEdges << endl;
  cgsize_t dim[2] = {6, numEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("Edges");
  cg_gopath(file, "/Base/Zone1/Edges");
  cg_array_write("EdgesData", CGNS_ENUMV(Integer), 2, dim, edges.data());

  // Write boundary edges
  // {p1, p2, c1, num}
  int numBoundaryEdges = boundaryEdges.size() / 5;
  cout << "Number of boundary edges: " << numBoundaryEdges << endl;
  cgsize_t boundaryDim[2] = {5, numBoundaryEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("BoundaryEdges");
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  cg_array_write("BoundaryEdgesData", CGNS_ENUMV(Integer), 2, boundaryDim, boundaryEdges.data());

  cgsize_t bcDim[2] = {3, 4};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("BCs");
  cg_gopath(file, "/Base/Zone1/BCs");
  cg_array_write("BCs", CGNS_ENUMV(Integer), 2, bcDim, bc_data);

  cg_close(file);
}
