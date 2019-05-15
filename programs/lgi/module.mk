$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
# module_units_f := 
module_programs_cpp := get_lgi_expansion get_spncci_seed_blocks
module_programs_cpp += transform_seed_blocks 
module_programs_cpp += generate_truncated_lsu3shell_basis
module_programs_cpp += test_lsu3shell get_u3s_subspaces Brel_subspace_lister
## module_programs_cpp += get_spncci_dimensions
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
