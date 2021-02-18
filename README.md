# INS-Solver
A Discontinuous Galerkin FEM CFD application. Solves the Incompressible Navier-Stokes equations.

Dependencies:
- OP2
- CGNS
- VTK
- cuBLAS
- OpenBLAS
- PETSc

Directory structure:
- The 'src' directory contains the code for the Incompressible Navier-Stokes solver.
  - The 'constants' directory contains precomputed matrices exported from MATLAB.
  - The 'cuBLAS' directory contains functions that perform the required BLAS operations on OP2 dats.
  - The 'openBLAS' directory is the same as above but for OpenBLAS.
  - The 'kernels' directory contains all the OP2 kernels.
  - All other directories are created by the OP2 code generation.
- The 'tools' directory contains an unstructured VTK to CGNS mesh conversion code.

Build instructions:
```
mkdir build
cd build
cmake ..
make
```

CMakeLists.txt sets the locations of libraries to what I have on my system but you can change these by passing the following flags to the CMake command: `OP2_DIR` and `CGNS_DIR`.

Creating CGNS grid:
```
cd build/tools
cp ../../naca0012.vtk .
./vtk2CGNS
```

Running Airfoil:
```
cd build/src
cp ../tools/naca0012.cgns .
./ins_cuda -iter 10 -alpha 0.0
./ins_openmp -iter 10 -alpha 0.0
./ins_seq -iter 10 -alpha 0.0
```

You can then use Paraview to view the result stored in `naca0012.cgns`.
