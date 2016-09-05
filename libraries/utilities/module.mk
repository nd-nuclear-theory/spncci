$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := utilities multiplicity_tagged
module_units_cpp-h := parsing
# module_units_f := 
module_programs_cpp := 
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

$(eval $(end-module))
