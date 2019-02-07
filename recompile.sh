#!/bin/bash
ifort -O2 -c kinds.f90
ifort -O2 -mkl -c lapack.f90
ifort -O2 -c minf.f90
ifort -O2 -mkl -c numFort.f90 -lmkl_lapack95_lp64
