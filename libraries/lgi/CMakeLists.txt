cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(lgi_headers
    lgi.h
    recurrence_lgi.h
    lgi_solver.h
    # null_solver.h
    lgi_unit_tensors.h
    # dimensions.h
    # su3rme.h
    # lgi_gen.h
    )
set(
  lgi_sources lgi.cpp lgi_solver.cpp
  lgi_unit_tensors.cpp 
  # null_solver.cpp dimensions.cpp su3rme.cpp lgi_gen.cpp
  )

# define am library and add source files
add_library(lgi ${lgi_sources})
add_library(spncci::lgi ALIAS lgi)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  lgi
  PUBLIC am::am
         basis::basis
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::sp3rlib
         spncci::u3shell
         spncci::lsu3shell
         spncci::spncci_basis
         OpenMP::OpenMP_CXX
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH lgi_PARENT_DIR)
target_include_directories(
  lgi INTERFACE $<BUILD_INTERFACE:${lgi_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

###############################################################################
# define installed targets
# actual installation is performed by top-level CMakeLists.txt
###############################################################################
install(
  TARGETS lgi
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${lgi_headers} DESTINATION include/spncci/lgi)


###############################################################################
# tests
###############################################################################
add_custom_target(lgi_tests)

add_executable(lgi_test EXCLUDE_FROM_ALL lgi_test.cpp)
target_link_libraries(
  lgi_test PRIVATE 
  fmt::fmt
  am::am
  spncci::utilities
  spncci::sp3rlib
  spncci::u3shell
  spncci::lsu3shell
  spncci::lgi
)
add_dependencies(lgi_tests lgi_test)


add_executable(recurrence_indexing_test EXCLUDE_FROM_ALL recurrence_indexing_test.cpp)

target_link_libraries(
  recurrence_indexing_test PRIVATE
  fmt::fmt
  am::am
  spncci::lgi
  OpenMP::OpenMP_CXX
  cppitertools::cppitertools
  spncci::spncci_basis
  spncci::u3shell
)


add_executable(recurrence_indexing_test_simple EXCLUDE_FROM_ALL recurrence_indexing_test_simple.cpp)

target_link_libraries(
  recurrence_indexing_test_simple PRIVATE
  fmt::fmt
  am::am
  spncci::lgi
  OpenMP::OpenMP_CXX
  cppitertools::cppitertools
  spncci::spncci_basis
  spncci::u3shell
  spncci::u3ncsm
)


if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests lgi_tests)
endif()

