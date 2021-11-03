#! /bin/bash

rm -rf Build
mkdir Build

cwd=$(pwd)/src

cd Build

ncols=$(tput cols)
barwidth=$((ncols/5))

echo -n '|'
ifort -O2 -c $cwd/kinds.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='
done
ifort -O2 -qmkl -c $cwd/lapack.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='
done
ifort -O2 -c $cwd/minf.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='
done
ifort -O2 -c $cwd/quadpack.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='
done
ifort -O2 -qmkl -c $cwd/NumFort.f90 -lmkl_lapack95_lp64
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='
done
echo '|'
cd ..
