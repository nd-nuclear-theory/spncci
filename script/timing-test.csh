#~/bin/csh
set Nmax = $1
set ranks = $2
echo "mpiexec --n ${ranks} /afs/crc.nd.edu/user/m/mcaprio/projects/lsu3shell/programs/tools/SU3RME_MPI model_space_03_03_Nmax${Nmax}.dat model_space_03_03_Nmax${Nmax}.dat relative_operators.dat > ~/spncci-timing/timing-Nmax${Nmax}-n${ranks}.txt"
mpiexec --n ${ranks} /afs/crc.nd.edu/user/m/mcaprio/projects/lsu3shell/programs/tools/SU3RME_MPI model_space_03_03_Nmax${Nmax}.dat model_space_03_03_Nmax${Nmax}.dat relative_operators.dat > ~/spncci-timing/timing-Nmax${Nmax}-n${ranks}.txt
