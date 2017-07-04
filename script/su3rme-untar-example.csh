  foreach run {mac0408 mac0411}
  cd /project/projectdirs/m2032/data/spncci/su3rme/run${run}
  foreach f (*.tgz)
    ${HOME}/code/spncci/script/su3rme-untar.csh $f ${SCRATCH}/data/spncci/su3rme/runmac0408
  end
  end
