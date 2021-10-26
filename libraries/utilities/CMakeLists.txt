cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(utilities_headers utilities.h nuclide.h)
#set(utilities_headers nuclide.h)
# set(utilities_sources "")  # no sources -- INTERFACE library

# define am library and add source files
add_library(utilities INTERFACE)
add_library(spncci::utilities ALIAS utilities)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(utilities INTERFACE basis::basis Eigen3::Eigen GSL::gsl am::am)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH utilities_PARENT_DIR)
target_include_directories(
  utilities INTERFACE $<BUILD_INTERFACE:${utilities_PARENT_DIR}>
                      $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS utilities
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${utilities_headers} DESTINATION include/spncci/utilities)
