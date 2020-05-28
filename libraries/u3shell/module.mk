
$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := tensor_labels relative_operator two_body_operator u3st_scheme
module_units_cpp-h += u3spn_scheme two_body_branching relative_branching
module_units_cpp-h += upcoupling unu3 unit_tensor_expansion
module_units_cpp-h += unit_tensor_space_u3s interaction_truncation
# module_units_f := 
module_programs_cpp += generate_unit_tensors_tb
module_programs_cpp_test := tensor_labels_test relative_operator_test two_body_operator_test 
module_programs_cpp_test += u3st_scheme_test u3spn_scheme_test upcoupling_test unu3_test
module_programs_cpp_test += unit_tensor_space_u3s_test relative_branching_test interaction_truncation_test
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
