################################################################
# directory trees
################################################################

# search prefix
#
#   additional path to search for required include and lib 
#   directories

search_prefix := $(GSL_DIR) $(BOOST_ROOT)

# directories to directly include in include path or lib path
#
#   only for use as a fallback if the traditional search prefix scheme
#   above fails for a given installation
search_dirs_include :=  $(EIGEN3_DIR)/include/eigen3
search_dirs_lib :=
## search_dirs_lib := /usr/local/gfortran/lib/x86_64 

# install prefix (only if you want to install the binaries somewhere)
#
# e.g., you would set to /usr/local to do a systemwide installation.
# This is analagous to the --prefix= option of autoconf installations.
install_prefix := install
################################################################
# C++ compiler-specific configuration
################################################################

# C++ compiler
CXX := mpicxx 
#CXX := icpc

# langage standard
CXXFLAGS += -std=c++14 -fopenmp -DCPP0X_STD_TR1

# avoid gcc 5 warnings on Eigen library
CXXFLAGS += -Wno-deprecated-declarations

# avoid gcc 5 warnings on shift operators
CXXFLAGS += -Wno-shift-count-overflow

# C++ compiler optimization and debugging
CXXFLAGS += -DNDEBUG -O3 #-Wall -W
#CXXFLAGS += -cxx=icpc -std=c++0x -shared-intel -DCPP0X_STD_TR1 -DMPICH_IGNORE_CXX_SEEK -openmp -O3
ifdef DEBUG
  CXXFLAGS += -g
endif

# parallel C++ compiler (DEPRECATED)
#   used in module.mk files as
#   program program.o: CXX := $(MPICXX)
##MPICXX := mpicxx -cxx=$(CXX)
##MPICXX := mpicxx
MPICXX := $(CXX)

################################################################
# FORTRAN compiler-specific configuration
################################################################

# FORTRAN compiler
# Example values:
#   for GCC 3.x: f77
#   for GCC 4.x: gfortran
#   for Intel: ifort
FC := mpif77
#FC := gfortran
#FC := ifort 
#FC := mpif77 -f77=ifort -recursive
#FC := mpif77 -f77=/sw/bin/gfortran
#FC := mpif77 -f77=gfortran -m64 -frecursive #frecursive needed for multithread safe su3lib
#FC := mpif90

# FORTRAN compiler flags
FFLAGS += -fopenmp -frecursive

# FORTRAN compiler optimization and debugging
FFLAGS += -O3
ifdef DEBUG
  FFLAGS += -g
endif

################################################################
# C++/FORTRAN linking 
#    with C++ main()
################################################################

# FORTRAN object libraries (added to LDLIBS)
# Example values, depending on the compiler you are using to compile
# the FORTRAN objects:
#   for GCC 3.x f77: -lg2c
#   for GCC 4.x gfortran: -lgfortran
#   for Intel ifort: -lifport -lifcore -limf
# A parallel fortran library might or might not be necessary, e.g.,
# not at NERSC, but with GCC 4.x and OpenMPI on the NDCRC:
#   -lmpi_mpifh -lgfortran

fortran_libs := -lgfortran
#fortran_libs += -lifport -lifcore -limf

# FORTRAN linking flags (added to LDFLAGS)
# Not yet needed but provided as hook.
fortran_flags := 

################################################################
# library configuration 
################################################################

# GNU OpenMP library
#
# needed for MFDn eigensolver with gcc 6

LDFLAGS += -lgomp

# SU3LIB numerical precision
#   Set flag SU3DBL for double precision or SU3QUAD for quad precision.
#   Note: quad precision requires ifort compiler

FFLAGS += -DSU3DBL
##FFLAGS += -DSU3QUAD
##FFLAGS += -DSU3QUAD_GNU
# machine-specific numerical library
# Gnu Scientific library
LDLIBS += -lgsl 
LDLIBS += -lgslcblas 
CPPFLAGS += -DHAVE_INLINE

# binary output using hdf5 
# LDLIBS += -lmascot -lz -lhdf5
#LDLIBS += -lz 

# boost
LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono

################################################################
# machine-dependent BLAS/LAPACK libraries (MKL/ACML)
#
# NOTE: This must be specified in a wrapper, not in this
# generic config.mk.
################################################################

# Intel Math Kernel library for MFDn eigensolver
# Needed when compiling with intel C++ & Intel Math Kernel library for MFDn eigensolver
#LDLIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -openmp

# ACML library for MFDn eigensolver
#LDLIBS += -lacml
#LDLIBS += -lirc
# ugly fix to force static linking of ACML
## LDLIBS += /afs/crc.nd.edu/x86_64_linux/scilib/acml/current/gfortran/gfortran64/lib/libacml.a

# generic LAPACK for easy system-independent compilation
#
# WARNING: The nonoptimized LAPACK is Not to be used for production!
## LDLIBS += /afs/crc.nd.edu/x86_64_linux/scilib/lapack/gnu/liblapack.a
