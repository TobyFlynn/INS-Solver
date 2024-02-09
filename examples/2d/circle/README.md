# Reinitalization of a circle

This is a test for the reinitalization algorithm. We have a perturbed level set
function exhibithing a wide range of gradients and the zero level set. 

The domain setup follows Ngu. [-2, 2], [2, 2] and the mesh refinement aims at
reconstructing the convergence figures for L1 norm.


## Running

Generate meshes
```bash
for n in $(seq 1 4)
do
    out_dir=m$n
    mkdir $out_dir
    scale=$((1.0/$n))
    gmsh \
        mesh.geo \
        -2 \
        -format vtk \
        -algo front2d \
        -nt 3 \
        -clscale $scale\
        -o $out_dir/mesh.vtk
    mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
        -file $out_dir/mesh.vtk \
        -bc none
    mv mesh.h5 $out_dir
done
```

```bash
ins2d_mpi -config config.ini -input m1/mesh.h5 -output m1/
```
