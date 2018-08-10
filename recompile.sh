#!/bin/bash

ifort -O2 -mkl -c kinds.f90
ifort -O2 -mkl -c pyplots.f90
ifort -O2 -mkl -c quadpack.f90
ifort -O2 -mkl -c lapack.f90
ifort -O2 -mkl -c numFort.f90
