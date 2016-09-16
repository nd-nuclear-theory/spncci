
$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := tensor_labels relative_operator two_body_operator u3st_scheme moshinsky
module_units_cpp-h += u3spn_scheme
module_units_cpp-h += import_interaction upcoupling unu3 unit_tensor_expansion

# module_units_f := 
module_programs_cpp := tensor_labels_test relative_operator_test two_body_operator_test 
module_programs_cpp += u3st_scheme_test moshinsky_test moshinsky_coefficient_test
module_programs_cpp += u3spn_scheme_test
module_programs_cpp += generate_unit_tensors_tb interaction_test upcoupling_test unu3_test #branching
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
