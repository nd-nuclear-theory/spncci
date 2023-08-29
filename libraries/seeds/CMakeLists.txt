cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(seeds_headers
  recurrence_seeds.h
    )
set(seeds_sources
  recurrence_seeds.cpp
  )

# define am library and add source files
add_library(seeds ${seeds_sources})
add_library(spncci::seeds ALIAS seeds)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  seeds
  PUBLIC am::am
         basis::basis
         fmt::fmt
         spncci::sp3rlib
         spncci::u3shell
         spncci::spncci_basis
         spncci::lgi
         spncci::u3ncsm
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH seeds_PARENT_DIR)
target_include_directories(
  seeds INTERFACE $<BUILD_INTERFACE:${seeds_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

###############################################################################
# define installed targets
# actual installation is performed by top-level CMakeLists.txt
###############################################################################
install(
  TARGETS seeds
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${seeds_headers} DESTINATION include/spncci/seeds)


###############################################################################
# tests
###############################################################################
add_custom_target(seeds_tests)


if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests seeds_tests)
endif()

add_executable(recurrence_seeds_test EXCLUDE_FROM_ALL recurrence_seeds_test.cpp)
target_link_libraries(
  recurrence_seeds_test PRIVATE
  fmt::fmt
  am::am
  spncci::utilities
  spncci::sp3rlib
  spncci::u3shell
  spncci::lsu3shell
  spncci::lgi
  spncci::u3ncsm
  spncci::seeds
)
