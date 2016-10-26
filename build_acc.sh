#! /bin/bash

# set up environment for compiling OpenACC program

module load CUDA/7.5.18 
module load PGI/16.3-GCC-4.9.3-2.25

export CXX=pgc++
export DEBUG=0

echo "CUDA loaded, PGI loaded, CXX set to pgc++ and DEBUG set to 0"
