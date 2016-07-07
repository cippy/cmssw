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

echo "Creating directory to store plots, if not yet existing ..."
mkdir -p plot
echo "Now launching executable ..."
echo "----------------------------"
echo " "
./EoverP $@
echo " "
echo "The end !"

