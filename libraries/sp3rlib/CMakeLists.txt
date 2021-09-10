cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(sp3rlib_headers multiplicity_tagged.h u3.h vcs.h sp3r.h u3coef.h
                    sp3r_operator.h
)
set(sp3rlib_sources u3.cpp vcs.cpp sp3r.cpp u3coef.cpp sp3r_operator.cpp)

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
         su3lib::su3lib
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
