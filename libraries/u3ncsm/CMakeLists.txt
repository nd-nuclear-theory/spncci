cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(u3ncsm_headers
    su3rme.h
    lgi_gen.h
    dimensions.h
    seed_gen.h
)
set(u3ncsm_sources
    su3rme.cpp
    lgi_gen.cpp
    dimensions.cpp
    seed_gen.cpp
)

# define spncci library and add source files
add_library(u3ncsm ${u3ncsm_sources})
add_library(spncci::u3ncsm ALIAS u3ncsm)

# ##############################################################################
# link dependencies
# ##############################################################################
target_link_libraries(
  u3ncsm
  PUBLIC am::am
         basis::basis
         fmt::fmt
         spncci::sp3rlib
         spncci::spncci_basis
         spncci::lgi
         OpenMP::OpenMP_CXX
         lsu3shell::SU3ME
         lsu3shell::LookUpContainers
         lsu3shell::LSU3
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH u3ncsm_PARENT_DIR)
target_include_directories(
  u3ncsm INTERFACE $<BUILD_INTERFACE:${u3ncsm_PARENT_DIR}>
                   $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS u3ncsm
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${u3ncsm_headers} DESTINATION include/spncci/u3ncsm)

# ##############################################################################
# define tests
# ##############################################################################


add_executable(dimensions_test EXCLUDE_FROM_ALL dimensions_test.cpp)
target_link_libraries(
  dimensions_test PRIVATE 
  fmt::fmt
  spncci::utilities
  spncci::lsu3shell
  spncci::u3ncsm
)

add_executable(modified_branching_test EXCLUDE_FROM_ALL modified_branching_test.cpp)
target_link_libraries(
  modified_branching_test PRIVATE
  fmt::fmt
  spncci::utilities
  spncci::lsu3shell
  spncci::u3ncsm
)

add_executable(print_lsu3shell_basis EXCLUDE_FROM_ALL print_lsu3shell_basis.cpp)
target_link_libraries(
  print_lsu3shell_basis PRIVATE
  fmt::fmt
  spncci::u3ncsm
)


# if(SPNCCI_MASTER_PROJECT)
#   add_dependencies(tests u3ncsm_tests)
# endif()
