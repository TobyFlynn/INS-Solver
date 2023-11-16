# 2D Example - Oscillating Droplet

## Generate Mesh
`gmsh mesh.geo -2 -format vtk -algo del2d -nt 2`
`mpirun -n 1 ../../../OP2-DG-Toolkit/bin/vtk2hdf5_2D -file mesh.vtk`

refRho = 1e-2
refVel = 1.0
refLen = 1.0
refMu  = 5e-5
refSurfTen = 0.5
mu0  = 1.0
mu1  = 100.0
rho0 = 1.0
rho1 = 100.0
