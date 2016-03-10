$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := u3coef
module_units_cpp-h := u3 vcs sp3r moshinsky

# module_units_f := 
module_programs_cpp := u3_test moshinsky_test u3coef_test
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
