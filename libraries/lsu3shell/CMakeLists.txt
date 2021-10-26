cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(lsu3shell_headers lsu3shell_operator.h lsu3shell_basis.h lsu3shell_rme.h su3rme.h)
set(lsu3shell_sources lsu3shell_operator.cpp lsu3shell_basis.cpp
                      lsu3shell_rme.cpp su3rme.cpp
)

# define am library and add source files
add_library(lsu3shell ${lsu3shell_sources})
add_library(spncci::lsu3shell ALIAS lsu3shell)

## find Tomas' lsu3shell libraries
find_package(lsu3shell REQUIRED)
# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  lsu3shell
  PUBLIC am::am
         basis::basis
         Boost::headers
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::utilities
         spncci::sp3rlib
         spncci::u3shell
         spncci::moshinsky
         OpenMP::OpenMP_CXX
         lsu3shell::libraries
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH lsu3shell_PARENT_DIR)
target_include_directories(
  lsu3shell INTERFACE $<BUILD_INTERFACE:${lsu3shell_PARENT_DIR}>
                      $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS lsu3shell
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${lsu3shell_headers} DESTINATION include/spncci/lsu3shell)

# ##############################################################################
# tests
# ##############################################################################

add_executable(lsu3shell_basis_test lsu3shell_basis_test.cpp)
target_link_libraries(
  lsu3shell_basis_test PRIVATE 
	lsu3shell::libraries
	fmt::fmt
	am::am
	spncci::utilities
	spncci::sp3rlib
	spncci::u3shell
	spncci::lsu3shell
)



