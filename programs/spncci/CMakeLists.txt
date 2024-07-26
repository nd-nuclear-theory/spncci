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
add_executable(kernels1 kernels1.cpp)
add_executable(kernels2 kernels2.cpp)
add_executable(kernels3 kernels3.cpp)
add_executable(TBDRME TBDRME.cpp)
add_executable(Sp3R_basis Sp3R_basis.cpp)
add_executable(dimension_calculator dimension_calculator.cpp)

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

target_link_libraries(
  kernels1
  PRIVATE spncci::spncci
)

target_link_libraries(
  kernels2
  PRIVATE spncci::spncci
)

target_link_libraries(
  kernels3
  PRIVATE spncci::spncci
)

target_link_libraries(
	TBDRME
  PRIVATE spncci::spncci
)

target_link_libraries(
        Sp3R_basis
  PRIVATE spncci::spncci
)

target_link_libraries(
        dimension_calculator
  PRIVATE spncci::spncci
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS spncci_exe kernels1 kernels2 kernels3 TBDRME Sp3R_basis dimension_calculator DESTINATION bin)
if(BUILD_LSU3SHELL)
  install(TARGETS generate_seeds_lsu3shell DESTINATION bin)
endif()
