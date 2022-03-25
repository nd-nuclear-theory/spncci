cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(u3shell_headers
    tensor_labels.h
    relative_operator.h
    two_body_operator.h
    u3st_scheme.h
    u3spn_scheme.h
    relative_branching.h
    upcoupling.h
    unu3.h
    unit_tensor_expansion.h
    unit_tensor_space_u3s.h
    interaction_truncation.h
)
set(u3shell_sources
    tensor_labels.cpp
    relative_operator.cpp
    two_body_operator.cpp
    u3st_scheme.cpp
    u3spn_scheme.cpp
    relative_branching.cpp
    upcoupling.cpp
    unu3.cpp
    unit_tensor_expansion.cpp
    unit_tensor_space_u3s.cpp
    interaction_truncation.cpp
)

# define am library and add source files
add_library(u3shell ${u3shell_sources})
add_library(spncci::u3shell ALIAS u3shell)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  u3shell
  PUBLIC am::am
         basis::basis
         Boost::headers
         mcutils::mcutils
         Eigen3::Eigen
         fmt::fmt
         spncci::sp3rlib
)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH u3shell_PARENT_DIR)
target_include_directories(
  u3shell INTERFACE $<BUILD_INTERFACE:${u3shell_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS u3shell
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${u3shell_headers} DESTINATION include/spncci/u3shell)


add_executable(relative_operator_test EXCLUDE_FROM_ALL relative_operator_test.cpp)
target_link_libraries(relative_operator_test PRIVATE
  spncci::spncci_basis
  fmt::fmt
  spncci::sp3rlib
  spncci::u3shell
)
