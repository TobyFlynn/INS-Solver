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
- The 'grids' directory contains some sample grids. The code is currently set up to deal with `cylinder.vtk` and `cylinder2.vtk`.

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
cp ../../grids/cylinder.vtk .
./vtk2CGNS -file cylinder.vtk -bc cylinder
```
The `-file` flag specifies the input file and the `-bc` flag specifies which boundary conditions to save with the CGNS file. Current valid options for `-bc` are `cylinder`, `poisson-test-0`, `poisson-test-1` and `vortex`.

Running INS:
```
cd build/src
cp ../tools/cylinder.cgns .
./ins_cuda -iter 1000 -pmethod 1 -input cylinder.cgns
./ins_openmp -iter 1000 -pmethod 1 -input cylinder.cgns
```

You can then use Paraview to view the end result which is stored in `end.cgns`. Flags for the INS solver are:
- `-iter` the number of iterations to run for
- `-pmethod` the method for Poisson solver (`0` for explicit matrix, `1` for matrix free)
- `-input` the input grid file
- `-output` the output directory (by default the same directory as the executable)
- `-save` save the flow every `x` iterations (currently only on the single node solvers)
- `-problem` the problem that the solver is solving (`0` for the cylinder problem, `1` for the vortex problem)
- `-multiphase` whether this is a one fluid (set flag to `0`, the default) or two fluid (set flag to `1`) problem

To run the vortex error checker:
```
cd build/src
./ins_cuda -iter 100 -pmethod 1 -input cylinder.cgns -problem 1
cd ../tools
cp ../src/end.cgns
./vortexError
```
The absolute error values at each grid point will then be saved to `err.cgns`.
