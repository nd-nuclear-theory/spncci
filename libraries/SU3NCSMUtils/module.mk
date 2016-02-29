$(eval $(begin-module))

################################################################
# unit definitions
################################################################

## module_units_h := CTuple combination 
module_units_h := CTuple
## module_units_cpp-h := CNullSpaceSolver clebschGordan factorial CRunParameters
# module_units_f := 
# module_programs_cpp :=

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

$(current-directory)/clebschGordan.o : $(current-directory)/factorial.h

$(eval $(end-module))
