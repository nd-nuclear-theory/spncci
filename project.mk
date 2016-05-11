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

modules := libraries/spncci libraries/u3shell libraries/sp3rlib libraries/utilities libraries/su3lib

# additional libraries -- imported
modules += libraries/cppformat
 
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
CPPFLAGS += -DHAVE_INLINE

# verbosity level
CPPFLAGS += -DNOVERBOSE -DNOVERBOSE_OMP

# program algorithm choices
#   C++ or FORTRAN wru3optimized
CPPFLAGS += -DNO_USE_LSU_WRU3 
#   hash function
CPPFLAGS += -DBOOSTHASH
#   map vs. hash unit tensor sectors 
CPPFLAGS += -DNOHASH_UNIT_TENSOR
#   precalculation and caching of U coefficients
CPPFLAGS += -DNOUSE_U_COEF_CACHE
#   map vs. hash for space lookup in BaseSpace
CPPFLAGS += -DINDEXING_BASE_HASH_SPACE

# debugging mode
CXXFLAGS += -g

#for lots of output
# -DVERBOSE 

# SU3LIB numerical precision
#   Set flag SU3DBL for double precision or SU3QUAD for quad precision.
#   Note: quad precision requires ifort compiler

FFLAGS += -DSU3DBL

# BOOST -- lsu3shell flags
## LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono

