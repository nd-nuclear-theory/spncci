################################################################
# project name
################################################################

project_name := spncci

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# libraries

# Caution: Order is important since used also in linking.  Caller must
# precede callee.

## modules := libraries/UNU3SU3 libraries/SU3NCSMUtils 

modules := libraries/spncci libraries/sp3rlib libraries/utilities libraries/su3lib

# additional libraries -- imported
## modules += libraries/cppformat  # need new gcc!
 
# additional libraries -- cloned as submodules
modules += libraries/am


#programs
modules += \
  programs/test 

################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras :=

################################################################
# additional project-specific make settings and rules
################################################################

# Gnu Scientific Library
LDLIBS += -lgsl 
LDLIBS += -lgslcblas 
CPPFLAGS += -DHAVE_INLINE -DBOOSTHASH -DHASH_UNIT_TENSOR -DVERBOSE

#for lots of output
# -DVERBOSE 

# to use hash table storage for unit tensor sectors
#-DHASH_UNIT_TENSOR

# SU3LIB numerical precision
#   Set flag SU3DBL for double precision or SU3QUAD for quad precision.
#   Note: quad precision requires ifort compiler

FFLAGS += -DSU3DBL

# BOOST
## LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono

