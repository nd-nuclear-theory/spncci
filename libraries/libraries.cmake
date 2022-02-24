cmake_minimum_required(VERSION 3.18)

# optionally use system-installed fmtlib
option(USE_SYSTEM_FMT "Use system-provided fmtlib" FALSE)

# ##############################################################################
# SU(3) coefficient library selection
# ##############################################################################

set(SPNCCI_SU3_LIBRARY_OPTIONS "su3wrc" "ndsu3lib" "SU3lib")
set(SPNCCI_SU3_LIBRARY "su3wrc" CACHE STRING "SU(3) coefficient library to use")
set_property(CACHE SPNCCI_SU3_LIBRARY PROPERTY STRINGS ${SPNCCI_SU3_LIBRARY_OPTIONS})

if(SPNCCI_SU3_LIBRARY STREQUAL "su3wrc")
  find_package(su3wrc REQUIRED)
  add_library(spncci::su3_library ALIAS su3wrc::su3wrc)
elseif(SPNCCI_SU3_LIBRARY STREQUAL "ndsu3lib")
  find_package(ndsu3lib REQUIRED)
  add_library(spncci::su3_library ALIAS ndsu3lib::ndsu3lib)
elseif(SPNCCI_SU3_LIBRARY STREQUAL "SU3lib")
  find_package(SU3lib REQUIRED)
  add_library(spncci::su3_library ALIAS SU3lib::SU3lib)
else()
  message(FATAL_ERROR "unknown SU(3) coefficient library ${SPNCCI_SU3_LIBRARY}")
endif()
message(STATUS "Found SU(3) coefficient library:  ${SPNCCI_SU3_LIBRARY}")


# ##############################################################################
# external dependencies
# ##############################################################################

find_package(Boost REQUIRED COMPONENTS headers)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Spectra REQUIRED)
find_package(su3lib REQUIRED)

## find Tomas' lsu3shell libraries
find_package(lsu3shell)

##find_package(GTest REQUIRED)
include(FetchContent)
FetchContent_Declare(
  GTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929 # v1.11.0
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(GTest)

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

# ##############################################################################
# define meta-library target "spncci::libraries"
# ##############################################################################

add_library(libraries INTERFACE)
add_library(spncci::libraries ALIAS libraries)
target_link_libraries(
  libraries
  INTERFACE spncci::utilities
            spncci::sp3rlib
            spncci::spncci_basis
            spncci::u3shell
            spncci::moshinsky
            spncci::lsu3shell
            spncci::lgi
            spncci::spncci
)

install(
  TARGETS libraries
  DESTINATION lib
  EXPORT spncciTargets
)
