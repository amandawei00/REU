#! /bin/bash

cd LHAPDF-6.2.1/

module load gcc/4.9.3
./configure --prefix=$PWD/..
cd ../bin/
export PATH=$PATH:$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
cd -

module load python/3.7.3
cd ../lib64/
export PATH=$PATH:$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
export PYTHONPATH=$PYTHONPATH:$PWD/python2.6/site-packages/
cd ../lib/
export PATH=$PATH:$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD

cd ..
