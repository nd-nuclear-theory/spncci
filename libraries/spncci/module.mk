$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := recurrence_indexing recurrence_indexing_spin
module_units_cpp-h := spncci_common spncci_basis
module_units_cpp-h += vcs_cache eigenproblem computation_control
module_units_cpp-h += explicit_construction recurrence parameters
module_units_cpp-h += io_control results_output decomposition transform_basis
module_units_cpp-h += hyperblocks_u3s variance
module_units_cpp-h += recurrence_indexing_spatial
#unit_tensor branching_u3lsj
# module_units_f :=
module_programs_cpp_test := spncci_basis_test
module_programs_cpp_test += recurrence_indexing_test
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
