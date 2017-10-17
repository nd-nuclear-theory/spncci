
$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h :=  moshinsky_xform relative_cm_xform
module_units_cpp-h += 

# module_units_f := 
module_programs_cpp :=  moshinsky_xform_test moshinsky_coefficient_test relative_to_twobody moshinsky_unitarity
module_programs_cpp +=  relative_cm_xform_test shell_comparison
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
