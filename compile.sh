#!/bin/bash

source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/setup.sh

cmake -B build
cmake --build build --target install --parallel
