#!/bin/csh

# Ex:
#   cd /project/projectdirs/m2032/data/spncci/su3rme/runmac0408
#   foreach f (*.tgz)
#     ${HOME}/code/spncci/script/su3rme-untar.csh $f $CSCRATCH/data/spncci/su3rme/runmac0408
#   end
set archive_name = $1
set target_parent_dir = $2

set basename = ${archive_name:r}
set target_dir = ${target_parent_dir}/${basename}
echo "${archive_name} -> ${target_dir}"

mkdir --parents ${target_dir}
nohup tar xf ${archive_name} --directory=${target_dir}
