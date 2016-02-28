################################################################
# project name
################################################################

project_name := spncci

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# libraries
# Caution: Order is important since used also in linking.
modules := libraries/UNU3SU3 libraries/SU3NCSMUtils libraries/su3lib

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
