$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h := spncci_common spncci_basis
module_units_cpp-h += branching_u3s branching_u3lsj
module_units_cpp-h += branching
module_units_cpp-h += unit_tensor vcs_cache
module_units_cpp-h += eigenproblem io_control parameters results_output decomposition
module_units_cpp-h += computation_control  # to split up
module_units_cpp-h += explicit_construction recurrence

# module_units_f := 
module_programs_cpp := spncci_basis_test branching_u3s_test
# module_programs_cpp += 
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
