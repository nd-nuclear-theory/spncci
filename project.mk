################################################################
# project name
################################################################

project_name := spncci

install_prefix := $(install_prefix)/spncci

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
modules += libraries/fmt

# additional libraries -- cloned as submodules
modules += libraries/basis libraries/am libraries/mcutils  # ordering note: mcutils is called by basis

#programs
modules += programs/operators programs/unit_tensors programs/su3calc programs/validation
#modules += programs/test
modules += programs/linear_algebra
modules += programs/spncci programs/lgi

modules += programs/rotor

################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras :=

################################################################
# project-specific make settings and rules
################################################################

# gsl
LDLIBS += -lgsl
LDLIBS += -lgslcblas
CPPFLAGS += -DHAVE_INLINE

# MKL and ScaLAPACK
CXXFLAGS += $(MKL_CXXFLAGS)
LDFLAGS += $(MKL_LDFLAGS)
LDLIBS += $(MKL_LDLIBS)
#LDLIBS += -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

# Eigen
CPPFLAGS += -DEIGEN_DONT_PARALLELIZE

# verbosity level
CPPFLAGS += -DNOVERBOSE -DNOVERBOSE_OMP

# basis submodule
#   map vs. hash for space lookup in basis module
CPPFLAGS += -DBASIS_HASH

# spncci program algorithm choices
#   map vs. hash unit tensor sectors
CPPFLAGS += -DNOHASH_UNIT_TENSOR
#   precalculation and caching of U coefficients
CPPFLAGS += -DUSE_U_COEF_CACHE

################################################################
# SU3LIB settings
################################################################

## # link to su3lib under lsu3shell
## LDFLAGS += -L../lsu3shell/libraries/su3lib

# numerical precision
#   SU3DBL: double precision
#   SU3QUAD: quad precision for ifort
#   SU3QUAD_GNU: quad precision for gnu gfortran
FFLAGS += -DSU3DBL

# compile-time setting of irrep range
#
#   SU3LM82: lambda+mu<82 instead of <42
FFLAGS += -DSU3LM82

################################################################
# spectra -- external template library
#
# This should arguably go in config.mk, but it is being placed
# here temporarily to keep config.mk standardized with the shell
# project.
################################################################

# % cd ${home}/code
# % git clone https://github.com/yixuan/spectra.git
# % setenv SPECTRA_DIR ${home}/code/spectra
#
# or
#
# % wget https://github.com/yixuan/spectra/archive/v0.5.0.tar.gz
#
# The spncci project uses the long form for eigen3 includes (e.g.,
# "eigen3/Eigen/Core"), but Spectra uses short form for eigen3
# includes (e.g., "Eigen/Core").  We therefore explicitly include the
# preprocessor option "-I${EIGEN3_DIR}/include/eigen3" via
# search_dirs_include.

search_prefix += $(SPECTRA_DIR)
search_dirs_include += $(EIGEN3_DIR)/include/eigen3
