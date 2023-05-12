$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
module_units_f := blocks conmat dbsr delta dlut drr3 dtu3r3 dwr3 ijfact multhy multu3 qtu3r3 \
u3mult wu3r3w xewu3 xewu3s xwu3 xwu3r3 wu39lm wzu3_optimized wru3_optimized \
# module_programs_cpp :=
# module_programs_f := su3genbk
# module_generated := fort.4

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

# rule for generating data file
# $(current-dir)/fort.4 : $(current-dir)/su3genbk
#	cd $(dir $@); ./$(notdir $<)

$(eval $(end-module))
