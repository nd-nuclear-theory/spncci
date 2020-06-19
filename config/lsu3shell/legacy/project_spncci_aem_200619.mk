################################################################
# project name
################################################################

project_name := su3shell

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# libraries
# Caution: Order is important since used also in linking.
modules := libraries/LSU3 libraries/SU3ME libraries/LookUpContainers libraries/UNU3SU3 libraries/SU3NCSMUtils libraries/eigensolver_MFDn libraries/su3lib

#programs
modules += programs/tools programs/upstreams \
  programs/LSU3shell \
  tests/LSU3 \
  programs/tools \
  programs/upstreams \
  programs/downstreams \
  programs/su3ncsm_hdf5IO \
  programs/su3ncsm_textIO \
  #programs/Densities \
################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras := su3shell-install.txt devel

################################################################
# additional project-specific make settings and rules
################################################################
