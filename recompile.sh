#!/bin/bash

rm -rf Build
mkdir Build

cwd=$(pwd)/src

cd Build

ncols=$(tput cols)
barwidth=$((ncols/3))

echo -n '|'
gfortran -c $cwd/kinds.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='    
done
gfortran -c $cwd/minf.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='    
done
gfortran -c $cwd/NumFortLite.f90
for ((i=1;i<=$barwidth;i++))
do
    echo -n '='    
done
echo '|'
cd ..
