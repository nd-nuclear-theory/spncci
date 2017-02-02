$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h := spncci_basis spncci_branching_u3s spncci_branching_u3lsj
module_units_cpp-h += unit_tensor

# module_units_f := 
module_programs_cpp := explicit_basis_construction
module_programs_cpp += spncci_basis_test spncci_branching_u3s_test
module_programs_cpp += unit_tensor_test
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
