$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := indexing_u3st moshinsky

# module_units_f := 
module_programs_cpp := indexing_u3st_test moshinsky_test
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