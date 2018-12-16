$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
# module_units_f := 
module_programs_cpp := spncci_MPI spncci rme2txt make_6j_table_q3
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
