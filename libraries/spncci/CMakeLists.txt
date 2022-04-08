cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(spncci_headers
    recurrence_indexing.h
    recurrence_indexing_spin.h
    recurrence_spatial.h
    spncci_common.h
    spncci_basis.h
    spncci_seeds.h
    vcs_cache.h
    eigenproblem.h
    computation_control.h
    explicit_construction.h
    recurrence.h
    parameters.h
    io_control.h
    results_output.h
    decomposition.h
    transform_basis.h
    hyperblocks_u3s.h
    variance.h
    recurrence_indexing_spatial.h
    recurrence_seeds.h
)
set(spncci_sources
    spncci_common.cpp
    spncci_basis.cpp
    spncci_seeds.cpp
    vcs_cache.cpp
    eigenproblem.cpp
    computation_control.cpp
    explicit_construction.cpp
    recurrence.cpp
    parameters.cpp
    io_control.cpp
    results_output.cpp
    decomposition.cpp
    transform_basis.cpp
    hyperblocks_u3s.cpp
    variance.cpp
    recurrence_indexing_spatial.cpp
    recurrence_seeds.cpp
    recurrence_spatial.cpp
)

# define spncci library and add source files
add_library(spncci ${spncci_sources})
add_library(spncci::spncci ALIAS spncci)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  spncci
  PUBLIC am::am
         basis::basis
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::sp3rlib
         spncci::u3shell
         spncci::lgi
         Spectra::Spectra
         OpenMP::OpenMP_CXX
         MPI::MPI_CXX
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH spncci_PARENT_DIR)
target_include_directories(
  spncci INTERFACE $<BUILD_INTERFACE:${spncci_PARENT_DIR}>
                   $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS spncci
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${spncci_headers} DESTINATION include/spncci/spncci)

# ##############################################################################
# define tests
# ##############################################################################

set(spncci_UNITS_TEST
    recurrence_indexing_test
    recurrence_mapping_test
    recurrence_seeds_test
    spncci_seeds_test
)

add_custom_target(spncci_tests)
foreach(test_name IN LISTS spncci_UNITS_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
  target_link_libraries(${test_name} PRIVATE spncci::spncci)
  add_dependencies(spncci_tests ${test_name})
endforeach()

target_link_libraries(
  recurrence_indexing_test PRIVATE fmt::fmt am::am spncci::lgi
                                   cppitertools::cppitertools spncci::spncci
)
target_link_libraries(
  recurrence_mapping_test PRIVATE fmt::fmt am::am spncci::lgi
                                  cppitertools::cppitertools spncci::spncci
)
target_link_libraries(
  recurrence_seeds_test PRIVATE fmt::fmt am::am spncci::lgi
                                cppitertools::cppitertools spncci::spncci
)

target_link_libraries(
  spncci_seeds_test
  PRIVATE fmt::fmt
          am::am
          spncci::lgi
          spncci::lsu3shell
          cppitertools::cppitertools
          spncci::spncci
          basis::basis
)

if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests spncci_tests)
endif()
