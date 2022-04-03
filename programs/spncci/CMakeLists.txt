cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################
if(BUILD_LSU3SHELL)
  add_executable(generate_seeds_lsu3shell generate_seeds_lsu3shell.cpp)
endif()

add_executable(spncci_exe spncci.cpp)
set_target_properties(spncci_exe PROPERTIES OUTPUT_NAME "spncci")
# note: spncci is already the name of the target for libraries/spncci, so we
# have to pick a different target name and override the output executable name

# ##############################################################################
# link dependencies
# ##############################################################################
if(BUILD_LSU3SHELL)
  message(STATUS "Including program/spncci/generate_seeds_lsu3shell in the build.")
  target_link_libraries(
    generate_seeds_lsu3shell PRIVATE
    fmt::fmt
    am::am
    spncci::utilities
    spncci::spncci_basis
    spncci::lgi
    spncci::u3ncsm
    spncci::spncci
  )
endif()

target_link_libraries(
  spncci_exe
  PRIVATE fmt::fmt
          OpenMP::OpenMP_CXX
          Spectra::Spectra
          am::am
          mcutils::mcutils
          spncci::lgi
          spncci::spncci
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS spncci_exe  DESTINATION bin)
if(BUILD_LSU3SHELL)
  install(TARGETS generate_seeds_lsu3shell DESTINATION bin)
endif()
