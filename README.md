# PyReweighting
A toolkit of python scripts "PyReweighting" is provided to facilitate the reweighting of accelerated molecular dynamics (aMD) simulations. PyReweighting implements a list of commonly used reweighting methods, including (1) exponential average that reweights trajectory frames by the Boltzmann factor of the boost potential and then calculates the ensemble average for each bin, (2) Maclaurin series expansion that approximates the exponential Boltzmann factor, and (3) cumulant expansion that expresses the reweighting factor as summation of boost potential cumulants.

# Usage
* Update dir_codes in the .sh files to the folder where you save the .py scripts

* Example command lines for using the .sh scripts:

./reweight-1d.sh $Emax $cutoff $binx $data $T

./reweight-2d.sh $Emax $cutoff $binx $biny $data $T

./reweight-3d.sh $Emax $cutoff $binx $biny $binz $data $T

Make sure set temperate for $T, which may have been ignored in previous calculations.

# Tutorial
* Prepare input file "weights.dat" in the following format:
Column 1: dV in units of kbT; column 2: timestep; column 3: dV in units of kcal/mol

For AMBER14:
awk 'NR%1==0' amd.log | awk '{print ($8+$7)/(0.001987*300)" " $2 " " ($8+$7)}' > weights.dat

For AMBER12:
awk 'NR%1==0' amd.log | awk '{print ($8+$7)" " $3 " " ($8+$7)*(0.001987*300)}' > weights.dat

For NAMD simulation:
grep "ACCELERATED MD" namd.log | awk 'NR%1==0' | awk '{print $6/(0.001987*300)" " $4 " " $6 " "$8}' > weights.dat

**./run.sh**

# NOTES
(1) Maclaurin series "pmf-2D-Phi_Psi-reweight-MC-order10-disc6.png" is equivalent to cumulant expansion on the 1st order "pmf-2D-Phi_Psi-reweight-CE1.xvg"
(2) Check out cumulant expansion to the 2nd order "pmf-2D-Phi_Psi-reweight-CE2.png"; normally it gives the most accurate result!
(3) The above python scripts work for any kind of reaction coordinates, e.g., atom distance, RMSD, or Principal Component Analysis (PCA) modes. You just need to change the default parameters "-Xdim -180 180 -discX 6 -Ydim -180 180 -discY 6" to the dimension and bin size of the corresponding reaction coordinates for reweighting.
(4) For aMD, the current reweighting scheme using cumulant expansion to the 2nd order is limited to aMD simulations of small systems, e.g., proteins with 10 - 40 residues. For larger proteins with more than 100 residues, the energetic noise would be too high for accurate reweighting.
(5) GaMD has been developed to reduce the energetic noise in simulations of large systems as noted in (4). Accurate energetic reweighting can be obtained from GaMD simulations of large biomolecular systems (e.g., with millions of atoms) using cumulant expansion to the 2nd order, provided that standard deviation of the boost potential is <= ~ 9 kcal/mol and a total of roughly 0.5 - 1 million frames are saved in the GaMD simulations.

# Citation
Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation, 10(7): 2677-2689.

Miao Y, Bhattarai, A and Wang, J (2020) Ligand Gaussian accelerated molecular dynamics (LiGaMD): Characterization of ligand binding thermodynamics and kinetics, Journal of Chemical Theory and Computation, 16(9): 5526–5547.
