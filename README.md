# RRKM Branching Ratios
Code and guide towards predicting branching ratios of reactions using RRKM theory.

### Calculating Rate Constants of Elementary Reactions with RRKM Theory

Here, we make use of the RRKM calculations made available from (Zhu,1993) in the fortran executable `rrkm.f` provided in the folder. This code makes use of an input file, in which a general-use version is given as `input.rrkm`. An example for how to fill out the input file is given in the executable `exampleRateConstantCalculation.sh`, with the execution and output provided below.

```bash
./exampleRateConstantCalculation.sh > exampleRateConstantCalculation.out
```

<details><summary>exampleRateConstantCalculation.out</summary>
<p>

```text
comment line 1 - Frequencies must be in cm^-1, Temperatures must be in Kelvin   
comment line 2 - Energies are relative to the intermediate and must be in kcal/m

EMIN =  58.10  KCAL/MOLE  DELE =   0.50  KCAL/MOLE  NE =    5

    5  RATE CONSTANTS TO BE CALCULATED     REACTION PATH DEGENERCY = 1.0


                    ***************************
                    *  DATA FOR THE MOLECULE  *
                    ***************************

NUMBER OF OSCILLATORS     NUMBER OF INTERNAL ROTORS

         24                         0

FREQUENCIES     1/CM

      3144.0      3133.0      3125.0      3063.0      3060.0      3045.0      1537.0
      1529.0      1516.0      1442.0      1402.0      1361.0      1347.0      1246.0
      1179.0      1157.0      1088.0      1014.0       970.0       938.0       912.0
       842.0       739.0       719.0

MOMENTS OF INERTIA     AMU*ANGSTROM**2

    174.00    174.00     73.00


ROTATIONAL TEMPERATURE IS    0.0000       J =     0


               ***********************************
               *  DATA FOR THE TRANSITION STATE  *
               ***********************************

*****STRUCTURE OF THE TRANSITION STATE IS SPECIFIED*****

POTENTIAL ENERGY     NUMBER OF OSCILLATORS     NUMBER OF INTERNAL ROTORS

    2.70                      23                         0

FREQUENCIES     1/CM

      3148.0      3139.0      3134.0      3066.0      3043.0      3031.0      1521.0
      1517.0      1513.0      1439.0      1416.0      1372.0      1334.0      1243.0
      1188.0      1166.0      1113.0      1045.0       982.0       943.0       886.0
       860.0       796.0

MOMENTS OF INERTIA     AMU*ANGSTROM**2

    181.00    181.00     69.00


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

RATE CONSTANTS CALCULATED USING THESE OPTIONS:     KROT = 0  NANH = 0  NCRIT = 0  JDEN = -1  JSUM = -1



    E        J        DENSITY        SUM        RATE CONSTANT

 KCAL/MOLE            1/CM-1                      SECOND-1


  58.10       0      0.1074D+08    0.3015D+10     0.8412D+13
  58.60       0      0.1201D+08    0.3385D+10     0.8449D+13
  59.10       0      0.1342D+08    0.3798D+10     0.8486D+13
  59.60       0      0.1498D+08    0.4259D+10     0.8523D+13
  60.10       0      0.1672D+08    0.4773D+10     0.8560D+13


```
    
The example above demonstrated how the rate constant changed over five different total energies. Try tinkering with it to see how it changes with the relative energy of intermediate to transition state, frequencies, rotational energy, and moments of inertia.

</p>
</details>

### Calculate Rate Constants Over All Elementary Steps in a Reaction

Consider the bimolecular reaction of vinylacetylene with the methylidyne radical to produce a hydrogen atom:

<p align="center">
    CH + C<sub>4</sub>H<sub>4</sub> ⟶ C<sub>5</sub>H<sub>4</sub> + H
</p>

There are a number of possible elementary steps: reactants adding or inserting across bonds, intermediates isomerizing between one another, and bond(s) breaking to form products. This code handles isomerizations and reactant/product formation with barriers, meaning any elementary step accompanied with a transition state can have its rate constant calculated. To do so, the intermediates' and transition states' energies, geometries, and frequencies must be known, also known as the potential energy profile or surface (PES) of the system. An example PES is provided for the example C<sub>5</sub>H<sub>5</sub> system in the `C5H5` directory:

```bash
$ tree C5H5
C5H5
├── PES.dat
├── geometries-frequenciesANDzpe.txt
└── geometries/
```

For this code, information must be stored as follows: the energies in `PES.dat`, the frequencies and zero-point energies (ZPEs) in `geometries-frequenciesANDzpe.txt`, and a list of XYZ files with optimized geometries of each intermediate and transition state in `geometries/`. After specifying the total energy and rotational excitation of the system in the script `getRateConstantMatrix.sh`, it may be executed to produce the following ouptut like so:

```bash
./getRateConstantMatrix.sh > rateConstantMatrix.dat
```

<details><summary>rateConstantMatrix.dat</summary>
<p>

```text
#i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 i23 i40 i41 i42 p1 p2 p3 p4 p5 p6 p7 p8 p9
                       i2           i3           i4           i5           i6           i7           i8           i9          i10          i11          i12          i13          i14          i15          i16          i17          i18          i19          i20          i21          i22          i23          i40          i41          i42           p1           p2           p3           p4           p5           p6           p7           p8           p9 
          i2          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.1370e+07   3.2660e+06   6.4850e+04          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   9.4850e+07 
          i3          0.0          0.0   8.6000e+11          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.2610e+09   2.4690e+09          0.0 
          i4          0.0   1.2370e+09          0.0          0.0          0.0          0.0          0.0   1.4940e+12          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          i5          0.0          0.0          0.0          0.0   1.4060e+11          0.0   2.4450e+09   3.8550e+09          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   6.6100e+08          0.0          0.0          0.0          0.0          0.0          0.0   1.8790e+10   2.4110e+10          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          i6          0.0          0.0          0.0   5.8960e+12          0.0   2.3910e+13          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          i7          0.0          0.0          0.0          0.0   3.0520e+10          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.0210e+08   2.1870e+09          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          i8          0.0          0.0          0.0   4.5010e+09          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.4520e+10          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          i9          0.0          0.0   2.2030e+12   1.0520e+08          0.0          0.0          0.0          0.0   1.4920e+08          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.7070e+09          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i10          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.0330e+09          0.0        error          0.0          0.0   4.1480e+11          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.1830e+10          0.0   1.1920e+11          0.0          0.0          0.0          0.0          0.0          0.0 
         i11          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0        error          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.2790e+10   1.4580e+11          0.0          0.0          0.0          0.0          0.0 
         i12          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.8810e+10          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.3170e+10          0.0          0.0          0.0   1.3510e+11          0.0          0.0          0.0          0.0          0.0          0.0 
         i13          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   5.2250e+08          0.0          0.0          0.0          0.0   1.3930e+07          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.1610e+10          0.0          0.0   6.4610e+07          0.0          0.0          0.0 
         i14          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   3.3860e+13          0.0          0.0          0.0          0.0   7.3150e+12          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i15          0.0          0.0          0.0          0.0          0.0   4.6810e+09          0.0          0.0          0.0          0.0          0.0          0.0   2.9780e+11          0.0        error          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   3.9710e+11          0.0          0.0          0.0          0.0          0.0 
         i16          0.0          0.0          0.0          0.0          0.0   1.1260e+12          0.0          0.0          0.0          0.0          0.0          0.0          0.0        error          0.0          0.0   2.9390e+07          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i17          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.8800e+08          0.0          0.0          0.0          0.0          0.0        error          0.0          0.0          0.0          0.0   4.7530e+09          0.0          0.0          0.0          0.0          0.0          0.0          0.0   3.7630e+08          0.0          0.0          0.0 
         i18          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.7840e+07          0.0          0.0        error          0.0   3.1940e+08          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i19          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0        error        error          0.0   5.4930e+08          0.0          0.0          0.0          0.0          0.0          0.0   4.3090e+10          0.0          0.0          0.0          0.0   8.2670e+08          0.0          0.0          0.0 
         i20          0.0          0.0          0.0   3.2690e+08          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   3.4590e+09          0.0   1.2530e+10   2.3080e+10          0.0          0.0          0.0          0.0   1.4790e+11          0.0          0.0          0.0   7.4230e+11          0.0          0.0          0.0          0.0 
         i21          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   2.6790e+09          0.0   5.5230e+12          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i22          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   9.1510e+10          0.0          0.0   2.3710e+11          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i23          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   0.0000e+00          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i40   1.7250e+08          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.0420e+13          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i41   5.6920e+07          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0   1.4320e+13          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
         i42   9.8030e+05          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p1          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p2          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p3          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p4          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p5          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p6          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p7          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p8          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
          p9          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0          0.0 
```

</p>
</details>

Here, the code makes a distinction between intermediates and products by setting rate constants of intermediate formation from products to zero, emulating low-pressure conditions.

In a system like above with many different pathways, checking each rate is difficult. The matrix of rate constants output, where each cell (i,j) is the rate of forming molecule i from molecule j, can be visualized as follows:

```bash
./visualizeRateConstantMatrix.sh rateConstantMatrix.dat rateConstantMatrix.png "RRKM Rate Constants Calculated Assuming 20.3 kJ/mol Total Energy"
```
