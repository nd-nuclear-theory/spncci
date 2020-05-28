$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := multiplicity_tagged
module_units_cpp-h := u3 vcs sp3r u3coef sp3r_operator sp3rcoef

# module_units_f := 
module_programs_cpp_test := u3_test sp3r_test u3coef_test vcs_test
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
