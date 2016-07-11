#!/bin/bash

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

outputDirName="plot/test_2016/"   # name of directory where files are stored (can be a single name or a path, they will be created from current directory below)
############################################
# WARNING  --> outputDirName must end with /
########################################### 
wwwBase="/afs/cern.ch/user/m/mciprian/www/"
wwwdir="/afs/cern.ch/user/m/mciprian/www/EoverP/"$outputDirName

echo "Creating directory to store output files, if not yet existing ..."
echo "mkdir -p $outputDirName"
mkdir -p $outputDirName
echo "Creating directory to store plots on website, if not yet existing ..."
echo "mkdir -p $wwwdir"
mkdir -p $wwwdir

# now must launch script to make plots visible from website
currentPath="$PWD"
cd $wwwBase
./copyphp.sh   # you should already have this script to be able to see files in your website                                                                           
cd $currentPath

echo "Now launching executable ..."
echo "----------------------------"
echo " "
echo "./EoverP $@ -dn $outputDirName"
./EoverP $@ -dn $outputDirName
echo " "
echo "Copying plots from ./$plotdir to $wwwdir"
echo "cp -r ${outputDirName}*.png $wwwdir"
cp -r ${outputDirName}*.png $wwwdir
echo "cp -r ${outputDirName}*.pdf $wwwdir"
cp -r ${outputDirName}*.pdf $wwwdir
echo "The end !"

