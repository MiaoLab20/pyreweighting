#!/bin/bash

# Core reweighting algorithms
nlen=20 # number of different simulation lengths to be used for free energy calculation
nruns=3 # number of simulations runs
T=300 # temperature
python PyReweighting-GaMD-TI.py -input gamd-ti-energy.dat -method stats -T $T -nruns $nruns -nlen $nlen | tee -a PyReweighting-GaMD-TI.log
python PyReweighting-GaMD-TI.py -input gamd-ti-energy.dat -method Gaussian -T $T -nruns $nruns -nlen $nlen | tee -a PyReweighting-GaMD-TI.log
# python PyReweighting-GaMD-TI.py -input gamd-ti-energy.dat -method ExpAverage -T $T -nruns $nruns -nlen $nlen | tee -a PyReweighting-GaMD-TI.log [For very small boost potentials]
# python PyReweighting-GaMD-TI.py -input gamd-ti-energy.dat -method noweight -T $T -nruns $nruns -nlen $nlen | tee -a PyReweighting-GaMD-TI.log [for cMD-TI]


