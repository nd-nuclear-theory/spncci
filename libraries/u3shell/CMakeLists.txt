cmake_minimum_required(VERSION 3.20)


# define units
set(u3shell_UNITS_H
  # recurrence_upcoupling
  operator_parameters
  operator_indexing_spin
  # operator_indexing_sectors

)
set(u3shell_UNITS_H_CPP
  tensor_labels
  relative_operator
  two_body_operator
  u3st_scheme
  u3spn_scheme
  relative_branching
  upcoupling
  unu3
  unit_tensor_expansion
  unit_tensor_space_u3s
  interaction_truncation
  operator_indexing_spatial

  )

# add headers and sources

set(u3shell_headers ${u3shell_UNITS_H} ${u3shell_UNITS_H_CPP})
list(TRANSFORM u3shell_headers APPEND .h)
set(u3shell_sources ${u3shell_UNITS_H_CPP})
list(TRANSFORM u3shell_sources APPEND .cpp)

# define am library and add source files
add_library(u3shell ${u3shell_sources})
add_library(spncci::u3shell ALIAS u3shell)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  u3shell
  PUBLIC am::am
         basis::basis
         Boost::headers
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::sp3rlib
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH u3shell_PARENT_DIR)
target_include_directories(
  u3shell INTERFACE $<BUILD_INTERFACE:${u3shell_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS u3shell
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${u3shell_headers} DESTINATION include/spncci/u3shell)


add_executable(relative_operator_test EXCLUDE_FROM_ALL relative_operator_test.cpp)
target_link_libraries(relative_operator_test PRIVATE
  spncci::spncci_basis
  fmt::fmt
  spncci::sp3rlib
  spncci::u3shell
)

add_executable(recurrence_upcoupling_test EXCLUDE_FROM_ALL recurrence_upcoupling_test.cpp)
target_link_libraries(recurrence_upcoupling_test PRIVATE
  spncci::spncci_basis
  fmt::fmt
  spncci::sp3rlib
  spncci::u3shell
)


  # add_executable(operator_indexin_test EXCLUDE_FROM_ALL operator_indexin_test.cpp)
  # target_link_libraries(operator_indexin_test PRIVATE
  #   spncci::spncci_basis
  #   basis::basis
  #   fmt::fmt
  #   spncci::sp3rlib
  # )
