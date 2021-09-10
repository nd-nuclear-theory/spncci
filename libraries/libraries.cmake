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
find_package(Spectra REQUIRED)
find_package(su3lib REQUIRED)

find_package(MKL)
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

# ##############################################################################
# define meta-library target "spncci::libraries"
# ##############################################################################

add_library(libraries INTERFACE)
add_library(spncci::libraries ALIAS libraries)
target_link_libraries(
  libraries
  INTERFACE libraries/su3lib
            libraries/utilities
            libraries/sp3rlib
            libraries/u3shell
            libraries/moshinsky
            libraries/lsu3shell
            libraries/lgi
            libraries/spncci
)

install(
  TARGETS libraries
  DESTINATION lib
  EXPORT spncciTargets
)
