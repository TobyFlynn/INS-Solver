# INS-Solver
A Discontinuous Galerkin FEM CFD application. Solves the Incompressible Navier-Stokes equations.

Dependencies:
- OP2
- CGNS
- VTK
- cuBLAS
- OpenBLAS
- PETSc
- Armadillo

Directory structure:
- The 'src' directory contains the code for the Incompressible Navier-Stokes solver.
  - The 'cuBLAS' directory contains functions that perform the required BLAS operations on OP2 dats.
  - The 'openBLAS' directory is the same as above but for OpenBLAS.
  - The 'kernels' directory contains all the OP2 kernels.
- The 'tools' directory contains an unstructured VTK to CGNS mesh conversion code.
- The 'grids' directory contains some sample grids. The grid to use for multiphase runs is `cylinder3.vtk`.
- The 'gen' directory is where the preprocessing script copies the code to before the OP2 code generation is called. The contents of this directory is a copy of the 'src' directory but with some constant values filled in (the OP2 code gen doesn't handle compiler definitions in certain places, hence this preprocessing step).

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
cp ../../grids/cylinder3.vtk cylinder.vtk
./vtk2CGNS -file cylinder.vtk -bc cylinder
```
The `-file` flag specifies the input file and the `-bc` flag specifies which boundary conditions to save with the CGNS file. Current valid options for `-bc` are `cylinder`, `poisson-test-0`, `poisson-test-1` and `vortex`. Use `cylinder` for multiphase runs with the `cylinder3.vtk` grid.

Running INS:
```
cd build/gen
cp ../tools/cylinder.cgns .
./ins_cuda -iter 1 -input cylinder.cgns -mu0 1.0 -mu1 10.0 -rho0 1.0 -rho1 10.0 -pre 1
```

You can then use Paraview to view the end result which is stored in `end.cgns`. Flags for the INS solver are:
- `-iter` the number of iterations to run for
- `-input` the input grid file
- `-output` the output directory (by default the same directory as the executable)
- `-save` save the flow every `x` iterations to `sol.cgns`
- `-mu0` and `-mu1` to set the relative dynamic viscosities of the fluids
- `-rho0` and `-rho1` to set the relative densities of the fluids
- `-pre` for whether or not to use preconditioning for the linear solves (set it to `1`)

Program structure:
- `ins.cpp` constains the main function of the program. Initialises a `Solver` object and then calls it for the specified number of iterations.
- `solver.cpp` contains the implementation of the `Solver` class. Contains `advection`, `pressure` and `viscosity` functions. When the `Solver` class is initialised, it loads the mesh and uses it to initalise a `DGMesh` object (from the OP2-DG-Toolkit library). It also initialises `INSData`, `LS`, `PressureSolve` and `ViscositySolve` objects (see below for details on each of these).
- `ins_data.cpp` contains the implementation of the `INSData` class. This class allocates OP2 dats needed for this application (`DGMesh` allocates some OP2 dats that are needed for general DG applications, but the `INSData` allocates OP2 dats for data specific to this application e.g. density.). It also initialises the value for some of these OP2 dats.
- `ls.cpp` contains the implementation of the `LS` class. This has OP2 dats specific to the level set method and functions for advecting the level set function and for reinitialisation (though this needs the troubled element detector to be implemented so isn't called at the moment).
- `poisson.cpp` contains implementations for the `PoissonSolve`, `PressureSolve` and `ViscositySolve` classes. `PressureSolve` and `ViscositySolve` extend the `PoissonSolve` class to set the appropriate factor for each solve and the correct preconditioner (AMG for pressure and Block-Jacobi for viscosity). `PoissonSolve` functions:
  - `solve` function called to solve the Poisson problem after all the setup functions have been called. Adds the BCs to the RHS and then calls PETSc's KSP solver (using a conjugate gradient solver).
  - `calc_rhs` function used with PETSc's shell matrix to perform matrix-vector multiplication (matrix free solve).
  - `precond` similar to `calc_rhs`, function used with PETSc's shell preconditioner to perform Block-Jacobi matrix free preconditioning.
  - `set_op` this function calculates the LHS matrix of the Poisson problem based on the current factor that has been set (for `PressureSolve` this is 1/rho and for `ViscositySolve` this is just nu)
- `poisson_cpu.cpp` and `poisson_gpu.cu` contain code for switching between OP2 dats and PETSc Vec/Mat objects (different for GPU and CPU).
