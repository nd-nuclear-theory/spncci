################################################################
#
# lsu3shell-customizations.mk
#
#################################################################
#
#  Background
# 
#  The configuration needed to build lsu3shell is similar to the generic
#  shell/spncci project configuration, defined by the config files under
#  spncci/config/ndconfig.  The main differences are:
#
#    - lsu3shell expects a different location for the Eigen library's header
#      files, so we must add:
#
#        search_dirs_include :=  $(EIGEN3_DIR)/include/eigen3
#
#    - Even though linking to gsl is not specific to the compile-time
#      environment, and therefore does not belong in config.mk, it is currently
#      missing from the lsu3shell project.mk and relegated to config.mk:
#
#        # GNU Scientific Library
#        LDLIBS += -lgsl 
#        LDLIBS += -lgslcblas 
#        CXXFLAGS += -DHAVE_INLINE
#
#    - lsu3shell involves MPI parallelized FORTRAN.  This may require extra
#      cluster-specific libraries (but soon the shell/spncci projects may need
#      this, as well!).  In particular, at the ndcrc (in config-gnu-ndcrc.mk):
#
#        fortran_libs := -lmpi_mpifh -lgfortran
#
#      Since config-gnu.mk already invokes -lgfortran, I think this can/should
#      be simplified to:
#
#        fortran_libs += -lmpi_mpifh
#
#    - lsu3shell uses boost.  Thus, search_prefix must include the Boost root directory:
#
#        search_prefix += $(BOOST_ROOT)
#
#      And specific Boost libraries must be linked:
#
#        LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono
#
#    - lsu3shell includes the MFDn eigensolver.  This seems to require the GNU
#      OpenMP library (to be revisited?):
#
#        LDFLAGS += -lgomp
#
#    - The su3lib precision needs to be set:
#
#        FFLAGS += -DSU3DBL
#
################################################################
#
# Description
#
# This file should contain customizations relative to the generic ndconfig
# config-<compiler>.mk files.  It should only contain customizations which are
# *generic* to any compiler environment.  Customizations for a specific compiler
# environment (e.g., gnu vs. intel, or at a specific cluster) should be kept
# separate.  The general calling scheme is that this file should be included
# from a wrapper file for each compiler/cluster:
#
#   config/lsu3shell/config-lsu3shell-<compiler>-<cluster>.mk:
#     include config/ndconfig/config-<compiler>-<cluster>.mk
#     include config/lsu3shell/lsu3shell-customizations.mk
#
# Then there may be other customizations needed for specific compiler
# environments, in particular, to support MPI (since for now ndconfig does not
# provide MPI support, though it will likely need to do so in the near future).
# E.g., at the ndcrc, we need:
#
#  # MPI configuration at ND CRC
#  CC := mpicxx
#  FC := mpif77
#  fortran_libs += -lmpi_mpifh
#
################################################################
#
# Setup
#
# The makefile include statements assume you have set up a symlink to the spncci
# config directory:
#
#    % cd ${HOME}/code/lsu3shell
#    % ln -s ${HOME}/code/spncci/config
#
################################################################
#
# 11/03/18 (mac): Extract lsu3shell-specific directives from spncci/config/lsu3shell/config-gnu.mk.
#
################################################################

# language standard
#
# LSU3shell requires "-std=c++14 -DCPP0X_STD_TR1".  Although ndconfig provides
# c++14, the TR1 standard is needed for some legacy code in lsu3shell.
#
# The origin of the NDEBUG flag is uncertain, but gcc 

CXXFLAGS += -DCPP0X_STD_TR1
CXXFLAGS += -DNDEBUG

# Eigen
#
# Provide special "include" path, since lsu3shell #include directives use a
# different path convention for the Eigen header files.

search_dirs_include +=  $(EIGEN3_DIR)/include/eigen3

# Boost library

search_prefix += $(BOOST_ROOT)
LDLIBS += -lboost_mpi -lboost_serialization -lboost_system -lboost_chrono

# GNU OpenMP library
#
# needed for MFDn eigensolver with gcc 6

LDFLAGS += -lgomp

# SU3LIB numerical precision
#
#   SU3DBL: double precision
#   SU3QUAD: quad precision under ifort
#   SU3QUAD_GNU: quad precision under gnu

FFLAGS += -DSU3DBL
##FFLAGS += -DSU3QUAD
##FFLAGS += -DSU3QUAD_GNU

# GNU Scientific Library
LDLIBS += -lgsl 
LDLIBS += -lgslcblas 
CPPFLAGS += -DHAVE_INLINE

################################################################
# special targets
################################################################

# target to generate just codes needed for spncci

base_programs = programs/tools/SU3RME_MPI

base_executables = $(addsuffix $(binary_ext),$(base_programs))

.PHONY: base
base: $(base_programs)

.PHONY: install-base
install-base: base
	@echo Installing base to $(install_dir_bin)...
	install -D $(base_executables) --target-directory=$(install_dir_bin)
