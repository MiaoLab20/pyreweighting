#!/bin/bash
dir_codes=~/codes/python/amd-2d-reweight
Emax=$1
cutoff=$2
binx=$3
biny=$4
binz=$5
data=$6
T=$7

echo "Usage: reweight-3d.sh $Emax $cutoff $binx $biny $binz $data $T"

# tryp
rb=7.5 # 6
ru=7.5 # 12

# CD
# rb=7.5 # 6
# ru=7.5 # 12

if [ -f weights.dat ]; then
echo "python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amdweight_CE -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-c1-$data.xvg pmf-3D-c1-$data-reweight-discx$binx-discy$biny-discz$binz.xvg
mv -v pmf-c2-$data.xvg pmf-3D-c2-$data-reweight-discx$binx-discy$biny-discz$binz.xvg
mv -v pmf-c3-$data.xvg pmf-3D-c3-$data-reweight-discx$binx-discy$biny-discz$binz.xvg
# mv -v 3D_Free_energy_surface.png pmf-3D-$data-reweight-CE2-discx$binx-discy$biny-discz$binz.png

fi # weights.dat

if [ -f exist.dat ]; then
echo "exist.dat"
echo "python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job noweight" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job noweight | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-3D-$data-noweight-discx$binx-discy$biny-discz$binz.xvg
mv -v 3D_Free_energy_surface.png pmf-3D-$data-noweight-discx$binx-discy$biny-discz$binz.png

if [ -f weights.dat ]; then
echo "python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amdweight -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amdweight -weight weights.dat | tee -a reweight_variable.log
mv -v pmf-$data.xvg pmf-3D-$data-reweight-exp-discx$binx-discy$biny-discz$binz.xvg
mv -v 3D_Free_energy_surface.png pmf-3D-$data-reweight-exp-discx$binx-discy$biny-discz$binz.png

echo "python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job histo -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job histo -weight weights.dat | tee -a reweight_variable.log
mv -v histo-$data.xvg histo-3D-$data-discx$binx-discy$biny-discz$binz.dat.xvg

echo "python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amd_dV -weight weights.dat" | tee -a reweight_variable.log
python $dir_codes/PyReweighting-3D.py -input $data -T $T -Emax $Emax -cutoff $cutoff -discX $binx -discY $biny -discZ $binz -lig_dG True -rb $rb -ru $ru -job amd_dV -weight weights.dat | tee -a reweight_variable.log
mv -v dV-stat-3D-$data.xvg dV-stat-3D-$data-reweight-discx$binx-discy$biny-discz$binz.xvg
mv -v dV-anharm-3D-$data.xvg dV-anharm-3D-$data-reweight-discx$binx-discy$biny-discz$binz.xvg

fi # weights.dat

fi

