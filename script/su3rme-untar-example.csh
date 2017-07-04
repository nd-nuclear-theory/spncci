#!/bin/csh
#SBATCH -p shared  
#SBATCH -n 1 
#SBATCH -t 12:00:00
#SBATCH -J untar

foreach run (mac0408 mac0411)
cd /project/projectdirs/m2032/data/spncci/su3rme/run${run}
foreach f (*.tgz)
  ${HOME}/code/spncci/script/su3rme-untar.csh $f ${CSCRATCH}/data/spncci/su3rme-expanded/run${run}
end
end
