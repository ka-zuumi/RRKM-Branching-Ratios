#!/bin/bash

# Example usage:
# ./getRateConstantMatrix.sh > rateConstantMatrix.dat

# The output of this script is a square matrix of numbers
# with rate constants of the isomerization going from the
# molecule row I to molecule column J ...
# the total rates of disappearance (row I to column I)
# are listed as zero

tmpfile=tmp

# This is a specially-pre-prepared file
rrkminputfile=input.rrkm

# File with energies (assumed in kJ/mol) and frequencies
# (assumed in cm^-1) and geometries (assumed in angstrom)
energyfile=C5H5/PES.dat
frequencyfile=C5H5/geometries-frequenciesANDzpe.txt
geometrydir=C5H5/geometries

# The total amount of energy, relative to the reactants (kJ/mol)
# We can think of this as the amount of excitation (collision
# energy, rotation, and vibration) given to the reactants
totalEnergy=20.3

# The rotational temperature; the RRKM code uses conservation
# laws (e.g. angular momentum) so needs to know the (average)
# amount of rotational energy
rotT=0

# i = 3N - 6
# j = 3N - 7
i=24
j=23

################################################################################################

# Compile the heavy-lifting RRKM script and the script that
# calculates moments of inertia
# Check here to make sure your fortran version can access
# LAPACK and BLAS and change compilation as necessary

gfortran -std=legacy rrkm.f -o rrkm.o
gfortran calculateMomentOfInertiaFromXYZ.f90 -o calculateMomentOfInertiaFromXYZ.o -llapack -lblas

################################################################################################

# Enumerate over all molecules found in the energy file
# that are intermediates or products
molecules=$(grep -v '^#' $energyfile | awk '{print $1}' | grep '^[ip]' | sort --version-sort | xargs)
TSs=$(grep -v '^#' $energyfile | awk '{print $1}' | grep '^[ts]' | sort --version-sort | xargs)

# List these molecules as a guide
echo "#$molecules"

# Add a header
printf "%12s " ""
echo "$molecules" | xargs -n1 | while read molecule; do
  printf "%12s " "$molecule"
done
printf "\n"

# For each molecule, get some information about it
echo "$molecules" | xargs -n1 | while read moleculeSTART; do

  # Get its geometry
  intXYZfile=$(ls $geometrydir/$moleculeSTART.xyz)

  # Get its energy (with ZPE)
  intE=$(grep "^ *$moleculeSTART " $energyfile | awk '{print $2}')

  # Get its frequencies
  intFreq=$(grep "^ *$moleculeSTART " $frequencyfile | awk '{$NF=""; $1=""; print $0}' | xargs -n1 | awk '{printf "%s,", $1} END {printf "\n"}')

  # If something is missing, exit out now
  if [ -z "$intXYZfile" ] || [ -z "$intE" ] || [ -z "$intFreq" ]; then
    echo "Big ERROR occurred for molecule '$moleculeSTART' found in $energyfile"
    echo "Exiting!"
    exit
  fi

  # Using the geometry file, calculate the moment of inertia
  # tensor in the principle axes frame
  intI=$(./calculateMomentOfInertiaFromXYZ.o < $intXYZfile | xargs -n1 | awk '{printf "%d,", $1} END {printf "\n"}')

  # Now see every other molecule the current molecule can isomerize to
  printf "%12s " "$moleculeSTART"
  echo "$molecules" | xargs -n1 | while read moleculeEND; do

    # If we are starting as a product, don't check...
    # we do not allow products to reform intermediates
    if echo "$moleculeSTART" | grep -q 'p'; then
      printf "%12s " "0.0"
      continue
    fi

    # Try two possible ways the transition state may
    # be named
    tsNAME1="ts_${moleculeSTART}_${moleculeEND}"
    tsNAME2="ts_${moleculeEND}_${moleculeSTART}"

    # If that transition state exists, get some information about it
    if echo "$TSs" | grep -q "$tsNAME1 \|$tsNAME2 "; then

      # Get the geomtry, energy, and frequencies, just like before
      tsXYZfile=$((ls $geometrydir/$tsNAME1.xyz; ls $geometrydir/$tsNAME2.xyz) 2> /dev/null | head -1)
      tsE=$(grep "^ *$tsNAME1 \|^ *$tsNAME2 " $energyfile | awk '{print $2}')
      tsFreq=$(grep "^ *$tsNAME1 \|^ *$tsNAME2 " $frequencyfile | awk '{$NF=""; $1=""; $2=""; print $0}' | xargs -n1 | awk '{printf "%s,", $1} END {printf "\n"}')

      # If something is missing, leave it as an "error" that we can
      # just put down as "no isomerization" potentially later
      if [ -z "$tsXYZfile" ] || [ -z "$tsE" ] || [ -z "$tsFreq" ]; then
        printf "%12s " "error"
	continue
      fi
      tsI=$(./calculateMomentOfInertiaFromXYZ.o < $tsXYZfile | xargs -n1 | awk '{printf "%d,", $1} END {printf "\n"}')

      # Then, prepare the RRKM input file
      energy1=$(echo "scale=2; ($totalEnergy - ($intE))*0.23900574/1" | bc -l)
      energy2=$(echo "scale=2; ($tsE - ($intE))*0.23900574/1" | bc -l)

      sed -e "s/lowestTotalEnergy/$energy1/" -e "s/totalEnergyIncrement/0.5/" -e "s/numberOfEnergyIncrements/5/" -e "s/rotationalTemperature/$rotT/" -e "s/intermediateNumberOfNormalModes/$i/" -e "s/intermediateFrequencies,/$intFreq/" -e "s/intermediateMomentOfInertiaTensor,/$intI/" -e "s/transitionStateNumberOfNormalModes/$j/" -e "s/transitionStateFrequencies,/$tsFreq/" -e "s/transitionStateMomentOfInertiaTensor,/$tsI/" -e "s/transitionStateEnergyWithZPE/$energy2/" $rrkminputfile > $tmpfile

      # Execute the RRKM code and print out the rate constants
      rateconstants=$(./rrkm.o < $tmpfile | awk '/DENSITY.*RATE CONSTANT/ {getline;getline;getline;getline;getline; while (NF==5) {print $5; getline}}' | xargs)

      # For now, just keep the first rate constant; convert all
      # "D"s to "E"s in scientific notation
      printf "%12s " $(echo "$rateconstants" | sed 's/D/E/g' | awk '{printf "%.4e\n", 1.0*$1}')

    # If the transition state doesn't exist, this rate is zero
    else
      printf "%12s " "0.0"
    fi
  done
  printf "\n"
done
