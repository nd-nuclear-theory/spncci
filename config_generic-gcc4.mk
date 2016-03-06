################################################################
# directory trees
################################################################

# search prefix
#   additional path to search for required libraries
#   (e.g., su3lib and eigen)
search_prefix := $(HOME)/local
search_dirs_include := 
search_dirs_lib :=

# install prefix
install_prefix := $(HOME)/local
# Note: You should reset to /user/local to do a systemwide 
# installation.  This is analagous to the --prefix= option of 
# autoconf installations.


################################################################
# C++ compiler-specific configuration
################################################################

# C++ compiler
CXX := g++

# static linking
# to reduce dependence on run-time library configuration changes
# (e.g., need to load same modules as at compile time)
CXXFLAGS += -static

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
