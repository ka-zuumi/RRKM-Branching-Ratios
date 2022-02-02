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

</p>
</details>

### Calculate Rate Constants Over All Elementary Steps in a Reaction




