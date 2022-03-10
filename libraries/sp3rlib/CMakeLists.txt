cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(sp3rlib_headers multiplicity_tagged.h u3.h u3boson.h vcs.h sp3r.h u3coef.h
                    sp3r_operator.h
)
set(sp3rlib_sources u3.cpp u3boson.cpp vcs.cpp sp3r.cpp sp3r_operator.cpp u3coef.cpp
                    u3coef_${SPNCCI_SU3_LIBRARY}.cpp
)

# define am library and add source files
add_library(sp3rlib ${sp3rlib_sources})
add_library(spncci::sp3rlib ALIAS sp3rlib)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  sp3rlib
  PUBLIC am::am
         basis::basis
         Boost::headers
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         GSL::gsl
         spncci::su3_library
         spncci::utilities
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH sp3rlib_PARENT_DIR)
target_include_directories(
  sp3rlib INTERFACE $<BUILD_INTERFACE:${sp3rlib_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS sp3rlib
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${sp3rlib_headers} DESTINATION include/spncci/sp3rlib)

# ##############################################################################
# define tests
# ##############################################################################

set(sp3rlib_UNITS_TEST
    sp3r_test
    u3boson_test
    u3_test
    u3coef_test
    vcs_test
    su3coupling_gtest
)

add_custom_target(sp3rlib_tests)
foreach(test_name IN LISTS sp3rlib_UNITS_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
  target_link_libraries(${test_name} spncci::sp3rlib)
  add_dependencies(sp3rlib_tests ${test_name})
endforeach()

target_link_libraries(u3coef_test OpenMP::OpenMP_CXX)
target_link_libraries(su3coupling_gtest GTest::gtest_main)

if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests sp3rlib_tests)
endif()
