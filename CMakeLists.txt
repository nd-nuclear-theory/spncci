cmake_minimum_required(VERSION 3.15)

# Determine if basis is built as a subproject (using add_subdirectory) or if it
# is the master project.
if(NOT DEFINED SPNCCI_MASTER_PROJECT)
  set(SPNCCI_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(SPNCCI_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

project(spncci CXX Fortran)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# configurable options
option(SPNCCI_U_COEFF_CACHE "Enable precalculation and caching of U coefficients" ON)
option(SPNCCI_HASH_UNIT_TENSOR "Use hash tables for unit tensor sectors" OFF)
option(BUILD_LSU3SHELL "Directly link to lsu3shell for seed calculation" ON)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

# compiler options
add_compile_options(
  $<$<CONFIG:Release>:-O3>
  $<$<CONFIG:Debug>:-O0>
  $<$<CONFIG:Debug>:-Wall>
  $<$<CONFIG:Debug>:-Wsign-conversion>
  $<$<CONFIG:Debug>:-Wextra>
  $<$<CONFIG:Debug>:-Wno-reorder>
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:-ggdb>
)

add_compile_definitions(
  EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
  )

if(SPNCCI_U_COEFF_CACHE)
  add_compile_definitions(USE_U_COEF_CACHE)
endif()
if(NOT SPNCCI_HASH_UNIT_TENSOR)
  add_compile_definitions(NOHASH_UNIT_TENSOR)
endif()

# make building tests optional TODO

# git revision
find_package(Git REQUIRED QUIET)
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --long --tags --always --dirty
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE PACKAGE_GIT_VERSION
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions(VCS_REVISION=${PACKAGE_GIT_VERSION})

################################################################################
# tests
################################################################################

if(SPNCCI_MASTER_PROJECT)
  add_custom_target(tests)
endif()

################################################################################
# include library and program setup
################################################################################
# include FetchContent declarations
include(FetchContent)
if(SPNCCI_MASTER_PROJECT)
  include(spncci_shopping_list.cmake)
endif()
# include libraries
include(libraries/libraries.cmake)

# include programs
include(programs/programs.cmake)

# ##############################################################################
# installation rules
# ##############################################################################

install(
  EXPORT spncciTargets
  NAMESPACE spncci::
  FILE spncciTargets.cmake
  DESTINATION lib/cmake/spncci
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/spncciConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/spncci"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/spncciConfigVersion.cmake"
  VERSION "${spncci_VERSION_MAJOR}.${spncci_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/spncciConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/spncciConfigVersion.cmake
        DESTINATION lib/cmake/spncci
)

export(
  EXPORT spncciTargets
  NAMESPACE spncci::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/spncciTargets.cmake"
)
