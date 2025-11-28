#!/bin/bash

rm -r build
mkdir build
cd build

cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
#make -j 8  > >(tee build.log) 2> >(tee error.log >&2)
ninja > >(tee build.log) 


cd ..
mv build/*.log ./
