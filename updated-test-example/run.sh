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
python PyReweighting-1D.py -input Psi.dat -cutoff 10 -Xdim -180 180 -disc 6 -Emax 20 -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-Psi.dat.xvg pmf-Psi-reweight-CE1.xvg
mv -v pmf-c2-Psi.dat.xvg pmf-Psi-reweight-CE2.xvg
mv -v pmf-c3-Psi.dat.xvg pmf-Psi-reweight-CE3.xvg

# Reweighting using Maclaurin series expansion
python PyReweighting-1D.py -input Psi.dat -disc 6 -Emax 20 -job amdweight_MC -order 10 -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Psi.dat.xvg pmf-Psi-reweight-MC-order10.xvg

# Reweighting using exponential average
python PyReweighting-1D.py -input Psi.dat -disc 6 -Emax 20 -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Psi.dat.xvg pmf-Psi-reweight.xvg

# NOTE: Check out cumulant expansion to the 2nd order "pmf-Psi-reweight-CE2.xvg"; normally it gives the most accurate result!

####################################################
# 2D data
# Prepare input data file "Phi_Psi" in two columns
# ptraj can be used for AMBER simulations

# Reweighting using cumulant expansion 
python PyReweighting-2D.py -cutoff 10 -input Phi_Psi -Xdim -180 180 -discX 6 -Ydim -180 180 -discY 6 -Emax 20 -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE1.xvg
mv -v pmf-c2-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE2.xvg
mv -v pmf-c3-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-CE3.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight-CE2.png

# Reweighting using Maclaurin series expansion
python PyReweighting-2D.py -input Phi_Psi -Emax 100 -discX 6 -discY 6 -job amdweight_MC -order 10 -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight-MC-order10-disc6.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight-MC-order10-disc6.png

# Reweighting using exponential average
python PyReweighting-2D.py -input Phi_Psi -Emax 20 -discX 6 -discY 6 -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-Phi_Psi.xvg pmf-2D-Phi_Psi-reweight.xvg
mv -v 2D_Free_energy_surface.png pmf-2D-Phi_Psi-reweight.png

# NOTES: 
# 1) Maclaurin series "pmf-2D-Phi_Psi-reweight-MC-order10-disc6.png" is equivalent to cumulant expansion on the 1st order "pmf-2D-Phi_Psi-reweight-CE1.xvg"
# 2) Check out cumulant expansion to the 2nd order "pmf-2D-Phi_Psi-reweight-CE2.png"; normally it gives the most accurate result!

