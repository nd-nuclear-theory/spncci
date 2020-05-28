$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
##module_units_cpp-h := 
module_units_cpp-h :=lgi lgi_solver null_solver lgi_unit_tensors

# module_units_f := 
module_programs_cpp_test := lgi_test null_solver_test
## module_programs_cpp += get_lgi_expansion

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
