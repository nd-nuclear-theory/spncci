$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
# module_units_f := 
module_programs_cpp := generate_lsu3shell_model_space
module_programs_cpp += generate_lsu3shell_relative_operators generate_lsu3shell_interaction
module_programs_cpp += generate_lsu3shell_two_body_unit_tensors
module_programs_cpp += explicit
## module_programs_cpp += check_two_body_unit_tensors 
## module_programs_cpp += check_lsu3shell_rmes
module_programs_cpp += generate_krel_squared generate_unit_tensor_statistics
## module_programs_cpp += lgi_solver_test
# module_programs_cpp += compute_unit_tensor_rmes

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

# $(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

# rule for generating data file
# $(current-dir)/fort.4 : $(current-dir)/su3genbk
#	cd $(dir $@); ./$(notdir $<)

$(eval $(end-module))
