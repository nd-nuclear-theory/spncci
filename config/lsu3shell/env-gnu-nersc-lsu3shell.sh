
if [[ ($NERSC_HOST == "edison") && ($CRAY_CPU_TARGET == "sandybridge") ]]; then
   # 9/12/17 (mac): fix erroneous default architecture sandybridge by reloading craype-ivybridge
   module load craype-ivybridge
fi

# different module names for eigen 3.3.3 on edison and cori as of 9/12/17 (mac)
if [[ $NERSC_HOST == "edison" ]]; then
   module load eigen
   export EIGEN3_DIR=${EIGEN_DIR}
elif [[ $NERSC_HOST == "cori" ]]; then
   module load eigen3
fi
## export EIGEN3_DIR=/global/project/projectdirs/m2032/opt/eigen-3.2.10
export SPECTRA_DIR=/global/project/projectdirs/m2032/opt/spectra-0.5.0

module swap PrgEnv-intel PrgEnv-gnu
module load cray-libsci
module load craype-hugepages2M
module load boost
module load gsl
module load python/3.6-anaconda-4.4
module load cray-hdf5-parallel
