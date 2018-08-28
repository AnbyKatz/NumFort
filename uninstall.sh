#!/bin/bash

tail -n +2 "makefile" > "makefile.tmp" && mv "makefile.tmp" "makefile"
tail -n +2 "makefile" > "makefile.tmp" && mv "makefile.tmp" "makefile"
tail -n +2 "makefile" > "makefile.tmp" && mv "makefile.tmp" "makefile"

rm recompile.sh
rm makefile
rm pyplots.py
rm customPlot.py
rm -r Pyplots

rm *.o
rm *.mod
