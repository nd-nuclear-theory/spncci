$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# Note: If we are only using the header files, no object files, must
# disable makefile library flag, so makefile does not attempt to
# create library archive.

## module_units_h := CTuple combination 
module_units_h := CTuple
## module_units_cpp-h := CNullSpaceSolver clebschGordan factorial CRunParameters
# module_units_f := 
# module_programs_cpp :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

$(current-directory)/clebschGordan.o : $(current-directory)/factorial.h

$(eval $(end-module))
