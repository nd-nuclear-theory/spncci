$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
# module_units_f := 
module_programs_cpp := generate_relative_u3st_operators truncate_interaction
module_programs_cpp += test_operator_symmetries
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
