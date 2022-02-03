#!/bin/bash

# Example usage:
# ./exampleRateConstantCalculation.sh

tmpfile=tmp

# This is a specially-pre-prepared file
rrkminputfile=input.rrkm

################################################################################################

# Compile the heavy-lifting RRKM script

gfortran -std=legacy rrkm.f -o rrkm.o

################################################################################################

# Example input:

# Total energy (potential + kinetic) available to the intermediate
# This energy is relative to the intermediate's energy in kcal/mol
totalEnergy=58.1

# Rotational temperature of the intermediate (necessary to account
# for angular momentum stuff in the RRKM code) in Kelvin
rotT=0

# If N is the number of atoms, then
#   i = number of positive frequencies for the intermediate (3N-6)
#   j = number of positive frequencies for the TS (3N-7)
i=24
j=23

# Intermediate and TS frequencies in cm^-1
intFreq="3144,3133,3125,3063,3060,3045,1537,1529,1516,1442,1402,1361,1347,1246,1179,1157,1088,1014,970,938,912,842,739,719,676,457,373,317,241,69,"
tsFreq="3148,3139,3134,3066,3043,3031,1521,1517,1513,1439,1416,1372,1334,1243,1188,1166,1113,1045,982,943,886,860,796,689,591,479,305,218,216,"

# Intermediate and TS energies in kcal/mol
intE=0.0
tsE=2.7

# Intermediate and TS moments of inertia in amu A^2
intI="174,174,73,"
tsI="181,181,69,"

################################################################################################

# Prepare the input for the input file

energy1=$(echo "$totalEnergy - ($intE)" | bc -l)
energy2=$(echo "$tsE - ($intE)" | bc -l)

sed -e "s/lowestTotalEnergy/$energy1/" -e "s/totalEnergyIncrement/0.5/" -e "s/numberOfEnergyIncrements/5/" -e "s/rotationalTemperature/$rotT/" -e "s/intermediateNumberOfNormalModes/$i/" -e "s/intermediateFrequencies,/$intFreq/" -e "s/intermediateMomentOfInertiaTensor,/$intI/" -e "s/transitionStateNumberOfNormalModes/$j/" -e "s/transitionStateFrequencies,/$tsFreq/" -e "s/transitionStateMomentOfInertiaTensor,/$tsI/" -e "s/transitionStateEnergyWithZPE/$energy2/" $rrkminputfile > $tmpfile

# Feed the input into the executable

./rrkm.o < $tmpfile #| tac | sed '/DENSITY.*RATE CONSTANT/q' | tac

                                                                    
