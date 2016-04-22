$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := sp_basis coef_cache unit_tensor

# module_units_f := 
<<<<<<< HEAD
module_programs_cpp := sp_basis_test coef_cache_test unit_tensor_test
=======
module_programs_cpp := sp_basis_test coef_cache_test unit_tensor_test
>>>>>>> 7ad888e298f97f36ce86b59baf7fde413604574e
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
