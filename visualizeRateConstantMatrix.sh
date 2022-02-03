#!/bin/bash

#
# Example usage:
# ./visualizeRateConstantMatrix.sh rateConstantMatrix.dat rateConstantMatrix.png "RRKM Rate Constants Calculated Assuming 20.3 kJ/mol Total Energy"
#

tmpfile1=tmp1
tmpfile2=tmp2

# For visualizing the energy
energyfile=C5H5/PES.dat # This has ZPE included

############################################################################

if [ $# -gt "3" ] || [ $# -lt 2 ]; then
  echo ""
  echo "Error: Requires exactly two or three arguments:"
  echo "       1) Rate constant matrix file"
  echo "       2) Name of output image file (must be .PNG)"
  echo "May have an optional third argument:"
  echo "       3) Text to put on image (to identify the graph)"
  echo "STOPPING"
  exit
fi

if [ $# -eq "3" ]; then
  gnuplotline=$3
else
  gnuplotline=""
fi

inputfile=$1
imagename=$2

############################################################################

# Get some metadata

# Get the number of molecules
Nmolecules=$(sed 1,2d $inputfile | wc -l)

# Mark 'error' labels with a negative rate constant
sed 's/error/ -1.0/g' $inputfile > $tmpfile1

# Get lower and upper bounds for the rate constants
sed 1,2d $tmpfile1 | awk '{for (i=2;i<=NF;i++) {if ($i > 0) {printf "%f\n", $i}}}' | sort -g > $tmpfile2
maxK=$(tail -1 $tmpfile2)
minK=$(head -1 $tmpfile2)

############################################################################

# Prepare the files for plotting

# tmpfile1
# Get rid of all headers and text information for the rate constants
sed 1,2d $tmpfile1 | awk '{$1=""; print $0}' > $tmpfile2
mv $tmpfile2 $tmpfile1

# tmpfile2
# Get energies ready with text labels for molecules
sed 1,2d $inputfile | awk '{print $1}' | while read molecule; do

   # If we can't find the molecules, leave its energy as "0"
  (echo "$molecule 0";
   grep "^ *$molecule  *[-0-9.]* *$" $energyfile) | tail -1

# Mark intermediates as "0" and products as "1"
done | awk '/^ *i/ {print $0, 0} /^ *p/ {print $0, 1}' > $tmpfile2

# Get some lower and upper bounds for the energy
maxE=$(awk '{printf "%f\n", $2}' $tmpfile2 | sort -g | tail -1)
minE=$(awk '{printf "%f\n", $2}' $tmpfile2 | sort -g | head -1)

############################################################################

gnuplot <<- EOF
set terminal pngcairo size 1400,1400
set output '$imagename'

set multiplot
set noborder
set noxtics
set noytics
unset key
set label "$gnuplotline " at screen 0.02,0.98
set label "Number of Molecules: $Nmolecules" at screen 0.02,0.02
set label "Red indicates an error for that calculation" tc rgb "red" at screen 0.98,0.02 right
set label "Gray indicates no transition state found or close-to-zero rate constant" tc rgb "gray" at screen 0.98,0.01 right

set bmargin at screen 0.05
set lmargin at screen 0.10
set rmargin at screen 0.85
set tmargin at screen 0.80

maxK = log10($maxK)
minK = log10($minK)
dK = (maxK - minK) / 10.0

set cbrange [minK-3*dK:maxK]
set cbtics font ",12"
set format cb "10^{%g}"
set cbtics add ("Error" minK-2*dK)
set cblabel "Rate Constant (s^{-1}) from ROW to COLUMN" font ",16" offset 2,0
set palette defined (minK-3*dK "red", minK-1.5*dK "red", minK-1.5*dK "gray", minK-0.5*dK "gray", minK-0.5*dK "black", maxK "green")
#plot '$tmpfile1' u 1:2:3 matrix w image
plot '$tmpfile1' u 1:2:((\$3)>0?(log10(\$3)):1/0) matrix w image,\
     '$tmpfile1' u 1:2:((\$3)==0?(minK-dK):1/0) matrix w image,\
     '$tmpfile1' u 1:2:((\$3)<0?(minK-2*dK):1/0) matrix w image

set border 3 lw 2

set bmargin at screen 0.85
set lmargin at screen 0.10
set rmargin at screen 0.85
set tmargin at screen 0.95

minE = $minE
maxE = $maxE
dE = (maxE - minE) / 10.0

set xrange [1-0.5:$Nmolecules+0.5]
unset xlabel
set yrange [minE-dE:maxE+dE]
set ylabel "Relative Energy\n(kJ/mol)" font ",16" offset -1,0
set ytics nomirror
set grid y lw 2

plot "$tmpfile2" u ((\$0)+1):((\$3)==0?(\$2):1/0) w p pt 7 lc rgb "black",\
     "$tmpfile2" u ((\$0)+1):((\$3)==1?(\$2):1/0) w p pt 9 lc rgb "blue",\
     "$tmpfile2" u ((\$0)+1):2:1 w labels offset 0,char 1 tc rgb "black"
unset multiplot
EOF

