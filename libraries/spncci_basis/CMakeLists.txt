cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(spncci_basis_headers
    recurrence_indexing.h
    recurrence_indexing_spin.h
    recurrence_indexing_spatial.h
)
set(spncci_basis_sources
    recurrence_indexing_spatial.cpp
)

# define spncci library and add source files
add_library(spncci_basis ${spncci_basis_sources})
add_library(spncci::spncci_basis ALIAS spncci_basis)

# ##############################################################################
# link dependencies
# ##############################################################################
target_link_libraries(
  spncci_basis
  PUBLIC am::am
         basis::basis
         fmt::fmt
         spncci::sp3rlib
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH spncci_basis_PARENT_DIR)
target_include_directories(
  spncci_basis INTERFACE $<BUILD_INTERFACE:${spncci_basis_PARENT_DIR}>
                   $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS spncci_basis
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${spncci_basis_headers} DESTINATION include/spncci/spncci_basis)

# ##############################################################################
# define tests
# ##############################################################################

set(spncci_basis_UNITS_TEST
    # recurrence_indexing_test
)

add_custom_target(spncci_basis_tests)
foreach(test_name IN LISTS spncci_basis_UNITS_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
  target_link_libraries(${test_name} PRIVATE spncci::spncci_basis)
  add_dependencies(spncci_basis_tests ${test_name})
endforeach()

# target_link_libraries(
#   recurrence_indexing_test PRIVATE fmt::fmt am::am spncci::lgi
#                                    cppitertools::cppitertools spncci::spncci
# )


if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests spncci_basis_tests)
endif()
