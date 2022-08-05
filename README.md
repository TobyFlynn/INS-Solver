# INS-Solver
A Discontinuous Galerkin FEM CFD application. Solves the Incompressible Navier-Stokes equations.

Dependencies:
- OP2 (for unstructure mesh operations)
- CGNS (INS mesh I/O)
- PETSc (for linear solves)
- VTK (Gmsh outputs mesh in VTK, a tool in the repo then converts this to a CGNS file)
- Armadillo (required by OP2-DG-Toolkit)
- cuBLAS (required by OP2-DG-Toolkit)
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
./ins_cuda -iter 1000 -input cylinder.cgns
./ins_openmp -iter 1000 -input cylinder.cgns
```

You can then use Paraview to view the end result which is stored in `end.cgns`. Flags for the INS solver are:
- `-iter` the number of iterations to run for
- `-input` the input grid file
- `-output` the output directory (by default the same directory as the executable)
- `-save` save the flow every `x` iterations in a file called `sol.cgns`

To run the interpolation tool:
```
cd build/tools
./interpolateFromCGNS -points /path/to/file_with_points_to_interpolate_to.txt -grid /path/to/original_grid.cgns -sol /path/to/output_file_from_ins.cgns -out /path/to/file_to_save_interpolated_data.txt
```

The format of the `file_with_points_to_interpolate_to.txt` file is a CSV file without headers:
```
0.0,0.0
1.0,0.5
0.2,0.75
```
