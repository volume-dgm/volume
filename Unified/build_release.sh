#!/bin/bash

rm -r release
mkdir release
cd release

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
      -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
make -j 8  > >(tee build.log) 2> >(tee error.log >&2)

echo "Done"

cd ..
cp .bin/* bin/
mv release/*.log ./
