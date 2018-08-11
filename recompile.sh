#!/bin/bash

ifort -O2 -mkl -c lapack.f90
ifort -O2 -mkl -c numFort.f90
