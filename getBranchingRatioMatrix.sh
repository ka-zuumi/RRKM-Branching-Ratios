#!/bin/bash

# Example usage:
# ./getBranchingRatioMatrix.sh rateConstantMatrix.dat > branchingRatioMatrix.dat
#
# To view the branching ratio of 9 products:
# sed 1,2d branchingRatioMatrix.dat | awk '{printf "%12s ", $1; for (i=NF-8;i<=NF;i++) {printf "%12s ", $i}; printf "\n"}'

tmpfile=tmp

inputfile=$1
correctedinputfile=errorRemoved-rateConstantMatrix.dat

###############################################################

# Compile the code that will calculate these ratios

gfortran numericallyCalculateBranchingRatio.f90 -o numericallyCalculateBranchingRatio.o

###############################################################

# Get some metadata
Nmolecules=$(sed 1,2d $inputfile | wc -l)
molecules=$(sed 1,2d $inputfile | awk '{print $1}' | xargs)

# Remove any potentially bad inputs
sed 's/error/ 0.0 /g' $inputfile > $correctedinputfile

# Print the same header
head -2 $inputfile

# Produce the output
./numericallyCalculateBranchingRatio.o $Nmolecules < <(sed 1,2d $correctedinputfile | awk '{$1=""; print $0}') > $tmpfile

# Add a 'side header'
paste -d' ' <(sed 1,2d $inputfile | awk '{printf "%12s\n", $1}') $tmpfile
