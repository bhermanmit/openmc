#!/bin/sh

set -ev

# For every CI, we will at least run phdf5-omp-debug

# Build MPICH
wget -q http://www.mpich.org/static/downloads/3.1.3/mpich-3.1.3.tar.gz
tar -xzvf mpich-3.1.3.tar.gz >/dev/null 2>&1
cd mpich-3.1.3
./configure --prefix=$PWD/../mpich_install -q
make -j -s >/dev/null 2>&1
make -s install >/dev/null 2>&1
cd ..
rm -rf mpich-3.1.3

# Build PHDF5
wget -q http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz
tar -xzvf hdf5-1.8.14.tar.gz >/dev/null 2>&1
cd  hdf5-1.8.14
CC=../mpich_install/bin/mpicc FC=../mpich_install/bin/mpif90 ./configure --prefix=$PWD/../phdf5_install -q --enable-fortran --enable-fortran2003 --enable-parallel
make -j -s >/dev/null 2>&1
make -s install >/dev/null 2>&1
cd ..
rm -rf hdf5-1.8.14

# Build HDF5 and PETSc for rest of debug tests
if [ "${TRAVIS_PULL_REQUEST}" = "true" ]; then

  # Build HDF5
  tar -xzvf hdf5-1.8.14.tar.gz >/dev/null 2>&1
  cd  hdf5-1.8.14
  CC=gcc FC=gfortran ./configure --prefix=$PWD/../hdf5_install -q --enable-fortran --enable-fortran2003
  make -j -s >/dev/null 2>&1
  make -s install >/dev/null 2>&1
  cd ..
rm -rf hdf5-1.8.14

fi