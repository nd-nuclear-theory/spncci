cmake_minimum_required(VERSION 3.18)

# optionally use system-installed fmtlib
option(USE_SYSTEM_FMT "Use system-provided fmtlib" FALSE)
# ##############################################################################
# external dependencies
# ##############################################################################

find_package(Boost REQUIRED COMPONENTS headers)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Spectra REQUIRED)

# ##############################################################################
# SU(3) coefficient library selection
# ##############################################################################

set(SPNCCI_SU3_LIBRARY_OPTIONS "su3wrc" "ndsu3lib" "SU3lib")
set(SPNCCI_SU3_LIBRARY "su3wrc" CACHE STRING "SU(3) coefficient library to use")
set_property(CACHE SPNCCI_SU3_LIBRARY PROPERTY STRINGS ${SPNCCI_SU3_LIBRARY_OPTIONS})

if(SPNCCI_SU3_LIBRARY STREQUAL "su3wrc")
  find_package(su3wrc)
  if(TARGET su3wrc::su3wrc)
    add_library(spncci::su3_library ALIAS su3wrc::su3wrc)
  else()
    FetchContent_MakeAvailable(su3wrc)
    add_library(spncci::su3_library ALIAS su3wrc)
  endif()

elseif(SPNCCI_SU3_LIBRARY STREQUAL "ndsu3lib")
  find_package(ndsu3lib)
  if(TARGET ndsu3lib::ndsu3lib)
    add_library(spncci::su3_library ALIAS ndsu3lib::ndsu3lib)
  else()
    set(SO3COEF_LIBRARY "wigxjpf")
    message(STATUS "building ndsu3lib with ${SO3COEF_LIBRARY} support")
    FetchContent_MakeAvailable(ndsu3lib)
    add_library(spncci::su3_library ALIAS ndsu3lib)
  endif()

elseif(SPNCCI_SU3_LIBRARY STREQUAL "SU3lib")
  find_package(SU3lib)
  if(TARGET SU3lib::SU3lib)
    add_library(spncci::su3_library ALIAS SU3lib::SU3lib)
  else()
    FetchContent_MakeAvailable(SU3lib)
    add_library(spncci::su3_library ALIAS SU3lib)
  endif()

else()
  message(FATAL_ERROR "unknown SU(3) coefficient library ${SPNCCI_SU3_LIBRARY}")

endif()

message(STATUS "Using SU(3) coefficient library:  ${SPNCCI_SU3_LIBRARY}")

if(BUILD_LSU3SHELL)
  message(STATUS "building spncci with lsu3shell support")
  ## find Tomas' lsu3shell libraries
  find_package(lsu3shell)
  if(NOT TARGET lsu3shell::lsu3shell)
    FetchContent_MakeAvailable(lsu3shell)
  endif()

endif()

if(NOT TARGET GTest::gtest_main)
  FetchContent_MakeAvailable(GTest)
endif()

#find_package(MKL) # disabled 21/10/31, broken on Mac (cvc)
if(TARGET MKL::MKL)
  target_link_libraries(Eigen3::Eigen INTERFACE MKL::MKL)
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
else()
  find_package(BLAS)
  if(TARGET BLAS::BLAS)
    target_link_libraries(Eigen3::Eigen INTERFACE BLAS::BLAS)
    target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_BLAS)
  endif()
endif()

option(EIGEN_DONT_PARALLELIZE "Turn off Eigen parallelization" ON)
if(EIGEN_DONT_PARALLELIZE)
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_DONT_PARALLELIZE)
endif()

# ##############################################################################
# internal dependencies
#
# note: since this script is include-d from the parent directory, we must
# prepend "library/" to each path
# ##############################################################################

if(USE_SYSTEM_FMT)
  find_package(fmt REQUIRED)
else()
  set(FMT_INSTALL ON)
  add_subdirectory(libraries/fmt)
endif()

add_subdirectory(libraries/am)
add_subdirectory(libraries/mcutils)
add_subdirectory(libraries/cppitertools)
add_subdirectory(libraries/basis)

# add_subdirectory(libraries/su3lib)
add_subdirectory(libraries/utilities)

add_subdirectory(libraries/sp3rlib)
add_subdirectory(libraries/u3shell)
add_subdirectory(libraries/moshinsky)
add_subdirectory(libraries/lsu3shell)
add_subdirectory(libraries/lgi)
add_subdirectory(libraries/spncci)
add_subdirectory(libraries/spncci_basis)

if(BUILD_LSU3SHELL)
  add_subdirectory(libraries/u3ncsm)
endif()

# ##############################################################################
# define meta-library target "spncci::libraries"
# ##############################################################################

add_library(spncci_libraries INTERFACE)
add_library(spncci::libraries ALIAS spncci_libraries)
target_link_libraries(
  spncci_libraries
  INTERFACE spncci::utilities
            spncci::sp3rlib
            spncci::spncci_basis
            spncci::u3shell
            spncci::u3ncsm
            spncci::moshinsky
            spncci::lsu3shell
            spncci::lgi
            spncci::spncci
)

If(BUILD_LSU3SHELL)
  target_link_libraries(spncci_libraries INTERFACE spncci::u3ncsm)
endif()

install(
  TARGETS spncci_libraries
  DESTINATION lib
  EXPORT spncciTargets
)
