#!/bin/bash

# Boot up
echo "Please make sure you have a working version of Python"
echo "with matplotlib and numpy plus ifort installed at"
echo "the bare minimum. Type C-c to exit, otherwise hit enter"
echo ""

read empty

# Remove files if they were previously installed
echo "Deleting files from previous installs"
rm -f pyplot.py
rm -f customPlot.py
echo "Compiling numFort library, may take a bit"
echo ""
./recompile.sh

####################################################################################

cwd=$(pwd)

touch pyplot.py
touch customPlot.py

cat .pyplotTemp >> pyplot.py
pypath=$(which python)
echo "#!$pypath" | cat - pyplot.py > temp && mv temp pyplot.py

cat .customTemp >> customPlot.py
echo "#!$pypath" | cat - customPlot.py > temp && mv temp customPlot.py

chmod +x pyplot.py
chmod +x customPlot.py

echo "NumFort wants to append somthing to your ~/.bashrc, make sure"
echo "the file exists"
read empty

if grep -Fxq "# Appended by NumFort, Path to directory" ~/.bashrc
then
    empty="nothing"
else
    echo "# Appended by NumFort, Path to directory" >> ~/.bashrc
    echo "export NumFortPath="$cwd/build/ >> ~/.bashrc
fi
    
source ~/.bashrc

echo ""
echo "Plotting code may be found in Pyplot directory"
echo "copy the created makefile template to your .f90 file."
echo "Type (use kinds and use numFort) in your .f90 file for base use"
echo "see readme documentation on specifics of available functions"
echo "Have fun! Let me know if they're any errors"
echo ""
