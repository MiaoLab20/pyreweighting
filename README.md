# PyReweighting
A toolkit of python scripts "PyReweighting" is provided to facilitate the reweighting of accelerated molecular dynamics (aMD) simulations. PyReweighting implements a list of commonly used reweighting methods, including (1) exponential average that reweights trajectory frames by the Boltzmann factor of the boost potential and then calculates the ensemble average for each bin, (2) Maclaurin series expansion that approximates the exponential Boltzmann factor, and (3) cumulant expansion that expresses the reweighting factor as summation of boost potential cumulants.

Usage:
* Update dir_codes in the .sh files to the folder where you save the .py scripts

* Example command lines for using the .sh scripts:
./reweight-1d.sh $Emax $cutoff $binx $data $T
./reweight-2d.sh $Emax $cutoff $binx $biny $data $T
./reweight-3d.sh $Emax $cutoff $binx $biny $binz $data $T
Make sure set temperate for $T, which may have been ignored in previous calculations.

Tutorial:
####################################################
# Prepare input file "weights.dat" in the following format: 
# Column 1: dV in units of kbT; column 2: timestep; column 3: dV in units of kcal/mol

# For AMBER14: 
# awk 'NR%1==0' amd.log | awk '{print ($8+$7)/(0.001987*300)" " $2 " " ($8+$7)}' > weights.dat

# For AMBER12: 
# awk 'NR%1==0' amd.log | awk '{print ($8+$7)" " $3 " " ($8+$7)*(0.001987*300)}' > weights.dat

# For NAMD simulation: 
# grep "ACCELERATED MD" namd.log | awk 'NR%1==0' | awk '{print $6/(0.001987*300)" " $4 " " $6 " "$8}' > weights.dat

####################################################
# 1D data
# Prepare input data file "Psi.dat" in one column, e.g., a dihedral angle Psi
# ptraj can be used for AMBER simulations

# Reweighting using cumulant expansion 
python PyReweighting-1D.py -input Psi.dat -T 300 -cutoff 10 -disc 6 -Emax 20 -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-Psi.dat.xvg pmf-Psi-reweight-CE1.xvg
mv -v pmf-c2-Psi.dat.xvg pmf-Psi-reweight-CE2.xvg
mv -v pmf-c3-Psi.dat.xvg pmf-Psi-reweight-CE3.xvg

# Reweighting using Maclaurin series expansion
python PyReweighting-1D.py -input Psi.dat -T 300 -disc 6 -Emax 20 -job amdweight_MC -order 10 -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Psi.dat.xvg pmf-Psi-reweight-MC-order10.xvg

# Reweighting using exponential average
python PyReweighting-1D.py -input Psi.dat -T 300 -disc 6 -Emax 20 -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Psi.dat.xvg pmf-Psi-reweight.xvg

# NOTE: Check out cumulant expansion to the 2nd order "pmf-Psi-reweight-CE2.xvg"; normally it gives the most accurate result!

####################################################
# 2D data
# Prepare input data file "Phi_Psi" in two columns
# ptraj can be used for AMBER simulations

# Reweighting using cumulant expansion 
python PyReweighting-2D.py -cutoff 10 -input Phi_Psi -T 300 -discX 6 -Ydim -180 180 -discY 6 -Emax 20 -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE1.xvg
mv -v pmf-c2-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE2.xvg
mv -v pmf-c3-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE3.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight-CE2.png

# Reweighting using Maclaurin series expansion
python PyReweighting-2D.py -input Phi_Psi -T 300 -Emax 100 -discX 6 -discY 6 -job amdweight_MC -order 10 -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-MC-order10-disc6.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight-MC-order10-disc6.png

# Reweighting using exponential average
python PyReweighting-2D.py -input Phi_Psi -T 300 -Emax 20 -discX 6 -discY 6 -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight.png

#################################################### # 
3D data 
# Prepare input data file "xyz" in 3 columns 
# cpptraj can be used for AMBER simulations 

# set "-lig_dG True" to calculate ligand binding free energy 
python PyReweighting-3D.py -input xyz -T 300 -Emax 100 -cutoff 500 -discX 1.0 -discY 1.0 -discZ 1.0 -lig_dG True -rb 7.5 -ru 7.5 -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log 
mv -v pmf-c1-xyz.xvg pmf-3D-c1-xyz-reweight-discx1.0-discy1.0-discz1.0.xvg 
mv -v pmf-c2-xyz.xvg pmf-3D-c2-xyz-reweight-discx1.0-discy1.0-discz1.0.xvg 
mv -v pmf-c3-xyz.xvg pmf-3D-c3-xyz-reweight-discx1.0-discy1.0-discz1.0.xvg  

NOTES: 
1) Maclaurin series "pmf-2D-Phi_Psi-reweight-MC-order10-disc6.png" is equivalent to cumulant expansion on the 1st order "pmf-2D-Phi_Psi-reweight-CE1.xvg"
2) Check out cumulant expansion to the 2nd order "pmf-2D-Phi_Psi-reweight-CE2.png"; normally it gives the most accurate result!
3) The above python scripts work for any kind of reaction coordinates, e.g., atom distance, RMSD, or Principal Component Analysis (PCA) modes. You just need to change the default parameters "-Xdim -180 180 -discX 6 -Ydim -180 180 -discY 6" to the dimension and bin size of the corresponding reaction coordinates for reweighting. 
4) The current reweighting scheme using cumulant expansion to the 2nd order is limited to aMD simulations of small systems, e.g., proteins with 10 - 40 residues. For larger proteins with more than 100 residues, the energetic noise would be too high for accurate reweighting. Further research has been focused on reducing such energetic noise in simulation of large proteins.

Citation:
Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation, 10(7): 2677-2689.

Miao Y, Bhattarai, A and Wang, J (2020) Ligand Gaussian accelerated molecular dynamics (LiGaMD): Characterization of ligand binding thermodynamics and kinetics, BioRxiv, https://doi.org/10.1101/2020.04.20.051979.
