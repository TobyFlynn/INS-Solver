# INS-Solver
A Discontinuous Galerkin FEM CFD application. Solves the Incompressible Navier-Stokes equations.

Dependencies:
- OP2
- PETSc
- HDF5
- INIPP
- VTK (required by OP2-DG-Toolkit)
- Armadillo (required by OP2-DG-Toolkit)
- OpenBLAS (required by OP2-DG-Toolkit)
- ParMETIS or PT-SCOTCH (required by OP2 and OP2-DG-Toolkit)

Directory structure:
- The 'src' directory contains the code for the Incompressible Navier-Stokes solver.
  - The 'kernels' directory contains all the OP2 kernels.
  - The 'io' directory contains code used to load the grid and save the solution.
  - The 'ls' directory contains code related to representing the interface between fluids using the level set method.
  - The 'poisson' directory contains code relating to solving the Poisson equation.
- The 'tools' directory contains an unstructured VTK to CGNS mesh conversion code.
- The 'grids' directory contains some sample grids.

Build instructions: see the example build script, `build.sh`, that demonstrates the preprocessing, OP2 translation, CMake and Make steps.

Running INS:
```
./ins_cuda -config path/to/config.ini -input path/to/grid.h5
```
