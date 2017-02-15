################################################################
# project name
################################################################

project_name := spncci

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# libraries

# Caution: Order of libraries is important since used also in linking.
# Calling library must precede callee library in listing.  That is,
# the most "basic" libraries must go last in the list (unless, of
# course, they are only template libraries, so nobody needs to link to
# them).

modules += libraries/spncci libraries/lgi libraries/lsu3shell 
modules += libraries/moshinsky libraries/u3shell libraries/sp3rlib  # ordering note: "mid-level" operations
modules += libraries/utilities libraries/su3lib  # ordering note: "low-level" operations, called by many other libraries

# additional libraries -- imported
modules += libraries/cppformat
 
# additional libraries -- cloned as submodules
modules += libraries/basis libraries/am libraries/mcutils  # ordering note: mcutils is called by basis

#programs
modules += programs/operators programs/unit_tensors programs/su3calc programs/validation
#modules += programs/test
modules += programs/spectra_test 

################################################################
# external template libraries
#
# This should arguably go in config.mk, but it is being placed
# here temporarily to keep config.mk standardized with the shell
# project.
################################################################

# % cd ${home}/projects
# % git clone https://github.com/yixuan/spectra.git
# % setenv SPECTRA_DIR ${home}/projects/spectra
#
# or 
#
# % wget https://github.com/yixuan/spectra/archive/v0.5.0.tar.gz
#
# The spncci project uses the long form for eigen3 includes (e.g.,
# "eigen3/Eigen/Core"), but Spectra uses short form for eigen3
# includes (e.g., "Eigen/Core").  We therefore explicitly include the
# preprocessor option "-I${EIGEN3_DIR}/include/eigen3".

search_prefix += $(SPECTRA_DIR)
search_dirs_include += $(EIGEN3_DIR)/include/eigen3


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
#   hash function
CPPFLAGS += -DBOOSTHASH
#   map vs. hash unit tensor sectors 
CPPFLAGS += -DNOHASH_UNIT_TENSOR
#   precalculation and caching of U coefficients
CPPFLAGS += -DUSE_U_COEF_CACHE
#   map vs. hash for space lookup in BaseSpace
CPPFLAGS += -DINDEXING_HASH

# debugging mode
#
# Define environment variable DEBUG on make command line to enable.
ifdef DEBUG
CXXFLAGS += -g
endif

# optimiation mode
CXXFLAGS += -O3

#for lots of output
# -DVERBOSE 

# SU3LIB numerical precision
#   Set flag SU3DBL for double precision or SU3QUAD for quad precision.
#   Note: quad precision requires ifort compiler

# double precision
# FFLAGS += -DSU3DBL
# quad precision for ifort
# FFLAGS += -DSU3QUAD
# quad precision for gnu gfortran
FFLAGS += -DSU3QUAD_GNU

# lambda+mu<82 instead of <42
FFLAGS += -DSU3LM82

# BOOST -- lsu3shell flags
## LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono

