cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(lgi_headers lgi.h lgi_solver.h null_solver.h lgi_unit_tensors.h)
set(lgi_sources lgi.cpp lgi_solver.cpp null_solver.cpp lgi_unit_tensors.cpp)

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

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS lgi
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${lgi_headers} DESTINATION include/spncci/lgi)
