cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(moshinsky_headers moshinsky_xform.h relative_cm_xform.h)
set(moshinsky_sources moshinsky_xform.cpp relative_cm_xform.cpp)

# define am library and add source files
add_library(moshinsky ${moshinsky_sources})
add_library(spncci::moshinsky ALIAS moshinsky)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  moshinsky
  PUBLIC am::am
         basis::basis
         Boost::headers
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::sp3rlib
         spncci::u3shell
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH moshinsky_PARENT_DIR)
target_include_directories(
  moshinsky INTERFACE $<BUILD_INTERFACE:${moshinsky_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS moshinsky
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${moshinsky_headers} DESTINATION include/spncci/moshinsky)
