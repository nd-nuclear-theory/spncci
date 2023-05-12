cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

# add_executable(explicit explicit.cpp)
add_executable(generate_krel_squared generate_krel_squared.cpp)
add_executable(
  generate_lsu3shell_interaction generate_lsu3shell_interaction.cpp
)
add_executable(
  generate_lsu3shell_model_space generate_lsu3shell_model_space.cpp
)
add_executable(
  generate_lsu3shell_relative_operators
  generate_lsu3shell_relative_operators.cpp
)
add_executable(
  generate_lsu3shell_two_body_unit_tensors
  generate_lsu3shell_two_body_unit_tensors.cpp
)
add_executable(
  generate_unit_tensor_statistics generate_unit_tensor_statistics.cpp
)
add_executable(generate_one_body_unit_tensor_labels generate_one_body_unit_tensor_labels.cpp) # added by J.H.
add_executable(generate_two_body_density_labels generate_two_body_density_labels.cpp) # added by J.H.

# ##############################################################################
# link dependencies
# ##############################################################################

# target_link_libraries(
#   explicit PRIVATE fmt::fmt mcutils::mcutils spncci::lgi spncci::spncci
# )
target_link_libraries(
  generate_krel_squared PRIVATE fmt::fmt spncci::lsu3shell spncci::sp3rlib
                                spncci::u3shell
)
target_link_libraries(
  generate_lsu3shell_interaction
  PRIVATE Boost::headers
          fmt::fmt
          OpenMP::OpenMP_CXX
          spncci::lsu3shell
          spncci::moshinsky
          spncci::sp3rlib
          spncci::u3shell
)
target_link_libraries(
  generate_lsu3shell_model_space PRIVATE fmt::fmt spncci::lsu3shell
                                         spncci::sp3rlib spncci::u3shell
)
target_link_libraries(
  generate_lsu3shell_relative_operators PRIVATE fmt::fmt spncci::lsu3shell
                                                spncci::sp3rlib spncci::u3shell
)
target_link_libraries(
  generate_lsu3shell_two_body_unit_tensors
  PRIVATE fmt::fmt spncci::lsu3shell spncci::sp3rlib spncci::u3shell
)
target_link_libraries(
  generate_unit_tensor_statistics PRIVATE fmt::fmt spncci::lgi
                                          spncci::lsu3shell spncci::u3shell
)
target_link_libraries(generate_one_body_unit_tensor_labels PRIVATE spncci::u3shell) # added by J.H.
target_link_libraries(generate_two_body_density_labels PRIVATE spncci::u3shell) # added by J.H.

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS generate_krel_squared
          generate_lsu3shell_interaction
          generate_lsu3shell_model_space
          generate_lsu3shell_relative_operators
          generate_lsu3shell_two_body_unit_tensors
          generate_unit_tensor_statistics
	  generate_one_body_unit_tensor_labels # added by J.H.
	  generate_two_body_density_labels # added by J.H.
          # explicit
  DESTINATION bin
)
