#!/bin/bash

# Boot up
echo "Please make sure you have a working version of Python"
echo "with matplotlib and numpy plus ifort installed at"
echo "the bare minimum. Type C-c to exit, otherwise hit enter"
echo ""

read empty

# Remove files if they were previously installed
echo "Deleting files from previous installs (May not exist)"
echo ""

rm recompile.sh
rm makefile
rm pyplots.py
rm customPlot.py
rm -r Pyplots

####################################################################################

touch recompile.sh
chmod +x recompile.sh

PLplotPath="/path/to/PLplotInstall"

echo ""
echo "did you install PLplot/Want to use shite plotting?"
echo "Please enter y or n"
read stringIn
while [ "$stringIn" != "y" ] &&  [  "$stringIn" != "n" ]
do
    echo "did you install PLplot?"
    echo "Please enter y or n"
    read stringIn
done
if [ "$stringIn" == "y" ]
then
    echo "Please input your path to PLplot in the form"
    echo "/path/to/PLplotInstallDirectory"
    read PLplotPath
fi

touch makefile
cat .makeTemp >> makefile

if [ "$stringIn" == "y" ]
then    
    rm makefile
    touch makefile
    cat .makeTempPL >> makefile
fi

cwd=$(pwd)
echo "DIR =" $cwd/ | cat - makefile > temp && mv temp makefile
echo "# DO NOT USE ~ FOR THE HOME DIRECTORY IN THIS PATH" | cat - makefile > temp && mv temp makefile
echo "# Directory to numFort" | cat - makefile > temp && mv temp makefile

echo ""

####################################################################################

echo "#!/bin/bash" >> recompile.sh
echo "ifort -O2 -c kinds.f90" >> recompile.sh
echo "ifort -O2 -mkl -c lapack.f90" >> recompile.sh

if [ "$stringIn" == "y" ]
then    
echo "ifort -O2 -c PLplots.f90 -I$PLplotPath/include/plplot -I$PLplotPath/lib/fortran/modules/plplot -L$PLplotPath/lib -lplplotfortran -lplplot" >> recompile.sh
fi

echo "ifort -O2 -mkl -c numFort.f90 -lmkl_lapack95_lp64" >> recompile.sh
echo "Compiling numFort library, may take a bit"
echo ""

./recompile.sh

####################################################################################

echo ""
echo "Please enter the version number for Python you have installed"
echo "Eg 3.5, 3.6, 3.7 are common versions"
read pyVersion

mkdir Pyplots
touch pyplots.py
touch customPlot.py
touch bashFortran.sh

cat .pyplotTemp >> pyplots.py
pypath=$(which python$pyVersion)
echo "#!$pypath" | cat - pyplots.py > temp && mv temp pyplots.py

cat .customTemp >> customPlot.py
echo "#!$pypath" | cat - customPlot.py > temp && mv temp customPlot.py

echo "#!/bin/bash" >> bashFortran.sh
echo "make all" >> bashFortran.sh
echo "./filename" >> bashFortran.sh
echo "python$pyVersion pyplots.py" >> bashFortran.sh
echo "# python$pyVersion customPlot.py" >> bashFortran.sh

chmod +x pyplots.py
chmod +x customPlot.py
chmod +x bashFortran.sh

cp pyplots.py Pyplots/pyplots.py
cp customPlot.py Pyplots/customPlot.py
cp makefile Pyplots/makefile
cp bashFortran.sh Pyplots/bashFortran.sh

echo ""
echo "Plotting code may be found in Pyplots directory"
echo "copy the created makefile template to your .f90 file."
echo "Type (use kinds and use numFort) in your .f90 file for base use"
echo "see readme documentation on specifics of available functions"
echo "Have fun! Let me know if they're any errors"
