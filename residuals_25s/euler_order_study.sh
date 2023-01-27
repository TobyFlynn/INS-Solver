#!/bin/bash

residual_dir=$(pwd)
number_iter=20000

cd ../build/gen

rm *.h5
./ins_cuda -input ../../grids/h5/vortex_square_1.h5 -iter $number_iter
mv residual.csv $residual_dir/residual_1.csv

rm *.h5
./ins_cuda -input ../../grids/h5/vortex_square_2.h5 -iter $number_iter
mv residual.csv $residual_dir/residual_2.csv

rm *.h5
./ins_cuda -input ../../grids/h5/vortex_square_3.h5 -iter $number_iter
mv residual.csv $residual_dir/residual_3.csv

rm *.h5
./ins_cuda -input ../../grids/h5/vortex_square_4.h5 -iter $number_iter
mv residual.csv $residual_dir/residual_4.csv
