# PyReweighting
A toolkit of python scripts "PyReweighting" is provided to facilitate the reweighting of accelerated molecular dynamics (aMD) simulations. PyReweighting implements a list of commonly used reweighting methods, including (1) exponential average that reweights trajectory frames by the Boltzmann factor of the boost potential and then calculates the ensemble average for each bin, (2) Maclaurin series expansion that approximates the exponential Boltzmann factor, and (3) cumulant expansion that expresses the reweighting factor as summation of boost potential cumulants.

# Usage
* Update dir_codes in the .sh files to the folder where you save the .py scripts

* Example command lines for using the .sh scripts:
./reweight-1d.sh $Emax $cutoff $binx $data $T

./reweight-2d.sh $Emax $cutoff $binx $biny $data $T

./reweight-3d.sh $Emax $cutoff $binx $biny $binz $data $T

Make sure set temperate for $T, which may have been ignored in previous calculations.

# Citation
Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation, 10(7): 2677-2689.

Miao Y, Bhattarai, A and Wang, J (2020) Ligand Gaussian accelerated molecular dynamics (LiGaMD): Characterization of ligand binding thermodynamics and kinetics, BioRxiv, https://doi.org/10.1101/2020.04.20.051979.
