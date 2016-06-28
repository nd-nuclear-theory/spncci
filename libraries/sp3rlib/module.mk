$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := u3 vcs sp3r u3coef lsushell_wru3 sp3r_operator un

# module_units_f := 
module_programs_cpp := u3_test sp3r_test u3coef_test vcs_test un_test
module_programs_cpp += su3_coupler  # to move to a utility program module
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
