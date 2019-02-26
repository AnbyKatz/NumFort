#!/bin/bash

# Boot up
echo "Please make sure you have a working version of Python"
echo "with matplotlib and numpy plus ifort installed at"
echo "the bare minimum. Type C-c to exit, otherwise hit enter"
echo ""

read empty

# Remove files if they were previously installed
echo "Deleting files from previous installs"
./uninstall.sh
echo "Compiling numFort library, may take a bit"
echo ""
./recompile.sh

####################################################################################

mkdir PythonPlot
cwd=$(pwd)
touch PythonPlot/Pyplot.py
touch PythonPlot/CustomPlot.py

cat Templates/PyplotTemp.py >> PythonPlot/Pyplot.py
pypath=$(which python)
echo "#!$pypath" | cat - PythonPlot/Pyplot.py > temp && mv temp PythonPlot/Pyplot.py

cat Templates/CustomTemp.py >> PythonPlot/CustomPlot.py
echo "#!$pypath" | cat - PythonPlot/CustomPlot.py > temp && mv temp PythonPlot/CustomPlot.py

chmod +x PythonPlot/Pyplot.py
chmod +x PythonPlot/CustomPlot.py

echo ""
echo "NumFort wants to append somthing to your ~/.bashrc, make sure"
echo "the file exists"
read empty

if grep -Fxq "# Appended by NumFort, Path to directory" ~/.bashrc
then
    empty="nothing"
else
    echo  >> ~/.bashrc
    echo "# Appended by NumFort, Path to directory" >> ~/.bashrc
    echo "export NumFortPath="$cwd/Build/ >> ~/.bashrc
        echo  >> ~/.bashrc
fi
    
source ~/.bashrc

echo ""
echo "Plotting code may be found in Pyplot directory"
echo "copy the created makefile template to your .f90 file."
echo "Type (use kinds and use numFort) in your .f90 file for base use"
echo "see readme documentation on specifics of available functions"
echo "Have fun! Let me know if they're any errors"
echo ""
