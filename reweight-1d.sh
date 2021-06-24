#!/bin/bash
dir_codes=~/codes/python/amd-2d-reweight
Emax=$1
cutoff=$2
binx=$3
data=$4
T=$5

echo "Usage: reweight-1d.sh $Emax $cutoff $binx $data $T"

if [ -f weights.dat ]; then
echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log"
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-$data.xvg pmf-c1-$data-reweight-disc$binx.dat.xvg
mv -v pmf-c2-$data.xvg pmf-c2-$data-reweight-disc$binx.dat.xvg
mv -v pmf-c3-$data.xvg pmf-c3-$data-reweight-disc$binx.dat.xvg

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-$data-noweight-disc$binx.dat.xvg

else

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-$data-noweight-disc$binx.dat.xvg

fi # weights.dat

if [ -f exist.dat ]; then
echo "exist.dat"

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amdweight -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-$data-amdweight-disc$binx.dat.xvg

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job histo -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job histo -weight weights.dat | tee -a reweight_variable.log
mv -v histo-$data.xvg histo-$data-disc$binx.dat.xvg

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amd_dV -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job amd_dV -weight weights.dat | tee -a reweight_variable.log

echo "python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-1D.py -input $data -T $T -disc $binx -Emax $Emax -cutoff $cutoff -job noweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-$data-noweight-disc$binx.dat.xvg

fi

