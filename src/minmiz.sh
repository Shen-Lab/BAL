# 1:current_directory 2:pdb_crd 3:charmm_setup 4:charmm_exe 5: crd_pdb


cd $1/minmiz/
cat $1/sampling/receptor.pdb $1/sampling/ligand.pdb > complex.pdb
$2 complex.pdb complex.crd

#/home/cyppsp/scripts/charmm_setup.prl complex.crd
$4  <minmiz.inp> charmm.out

$5 receptor.crd receptor.pdb
$5 ligand.crd ligand.pdb

