$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
# module_units_f := 
module_programs_cpp := triaxial_rotor sp3r_rotor
module_programs_cpp +=  transitions
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
