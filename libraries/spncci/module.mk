$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h := spncci_common spncci_basis 
module_units_cpp-h += branching_u3s branching_u3lsj branching branching2
module_units_cpp-h += unit_tensor vcs_cache eigenproblem computation_control
module_units_cpp-h += explicit_construction recurrence parameters
module_units_cpp-h += io_control results_output decomposition transform_basis
module_units_cpp-h += hyperblocks_u3s

# module_units_f := 
module_programs_cpp := spncci_basis_test
module_programs_cpp += branching2_test 
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
