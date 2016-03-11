################################################################
# directory trees
################################################################

# search prefix
#
#   additional path to search for required include and lib 
#   directories
#
# e.g., if the eigen include files are
#
#   $(HOME)/local/opt/eigen-3.0.3/include/eigen3/Eigen/{Array,Cholesky,...,Eigen,...}
#
# then the search prefix list should include
#
#   $(HOME)/local/opt/eigen-3.0.3/

search_prefix := $(HOME)/local/opt/eigen-3.0.3/

# directories to directly include in include path or lib path
#
#   only for use as a fallback if the traditional search prefix scheme
#   above fails for a given installation
search_dirs_include := 
search_dirs_lib :=

# install prefix (only if you want to install the binaries somewhere)
#
# e.g., you would set to /usr/local to do a systemwide installation.
# This is analagous to the --prefix= option of autoconf installations.
install_prefix := 


################################################################
# C++ compiler-specific configuration
################################################################

# C++ compiler
CXX := g++

# static linking
# to reduce dependence on run-time library configuration changes
# (e.g., need to load same modules as at compile time)
#
# CAVEAT: may cause trouble with system library linkage on OS X
CXXFLAGS += -static

# C++11 standard
CXXFLAGS += -std=c++0x

# parallel C++ compiler
#   used in module.mk files as
#   program program.o: CXX := $(MPICXX)
MPICXX := mpicxx -cxx=$(CXX)

################################################################
# FORTRAN compiler-specific configuration
################################################################

# FORTRAN compiler 
# Example values:
#   for GCC 3.x: f77
#   for GCC 4.x: gfortran
#   for Intel: ifort
FC := gfortran

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
fortran_libs := -lgfortran

# FORTRAN linking flags (added to LDFLAGS)
# Not yet needed but provided as hook.
fortran_flags :=
