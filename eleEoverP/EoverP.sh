#!/bin/bash

plotdir="plot/2016_Et"
rootdir="output/2016_Et"
wwwBase="/afs/cern.ch/user/m/mciprian/www/"
wwwdir="/afs/cern.ch/user/m/mciprian/www/EoverP/"$plotdir


echo
echo "Compiling EoverP.C ..."
g++ -Wall -pedantic -lm -o EoverP EoverP.C `rootlib`

if [ $? -ne 0 ]
then
    echo "Compilation failed!"
    return 0
fi

echo "Compilation succeeded :)"

for option in "$@";
do
    if [ $option = "-c" ]
    then
    # if option -c passed, we just want to complile and exit script
	return 0
    fi
done

echo "Creating directory to store plots, if not yet existing ..."
mkdir -p $plotdir
echo "Creating directory to store root files, if not yet existing ..."
mkdir -p $rootdir
echo "Creating directory to store plots on website, if not yet existing ..."
mkdir -p $wwwdir

# now must launch script to make plots visible from website
currentPath="$PWD"
cd $wwwBase
./copyphp.sh   # you should already have this script to be able to see files in your website                                                                           
cd $currentPath

echo "Now launching executable ..."
echo "----------------------------"
echo " "
./EoverP $@
echo " "
echo "Copying plots from ./$plotdir to $wwwdir"
cp -r ${plotdir}/*.png $wwwdir
cp -r ${plotdir}/*.pdf $wwwdir
echo "The end !"

