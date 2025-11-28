#!/bin/bash

cd build
#make -j 8  > >(tee build.log) 2> >(tee error.log >&2)
ninja > >(tee build.log)

cd ..
mv build/*.log ./
