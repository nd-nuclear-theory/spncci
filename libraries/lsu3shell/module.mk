$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
##module_units_cpp-h := lsu3shell_interface
module_units_cpp-h := lsu3shell_operator lsu3shell_basis lsu3shell_rme

# module_units_f := 
module_programs_cpp_test := lsu3shell_basis_test lsu3shell_rme_test

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
