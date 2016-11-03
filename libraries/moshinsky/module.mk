
$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h :=  moshinsky
module_units_cpp-h += 

# module_units_f := 
module_programs_cpp :=  moshinsky_test moshinsky_coefficient_test relative_to_twobody moshinsky_unitarity
module_programs_cpp +=  relative_cm_test
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
