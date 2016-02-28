$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := UNU3SU3Basics CUNMaster CSU3Master
# module_units_f := 
# module_programs_cpp :=

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

# $(current-dir)/UNU3SU3Basics.o: libraries/SU3NCSMUtils/CTuple.h
# $(current-dir)/CSU3Master.o: $(current-dir)/UNU3SU3Basics.h libraries/SU3NCSMUtils/CTuple.h

$(eval $(end-module))
