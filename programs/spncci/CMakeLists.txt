cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

add_executable(make_6j_table_q3 make_6j_table_q3.cpp)
add_executable(generate_seeds_lsu3shell generate_seeds_lsu3shell.cpp)
add_executable(spncci_exe spncci.cpp)
set_target_properties(spncci_exe PROPERTIES OUTPUT_NAME "spncci")
# note: spncci is already the name of the target for libraries/spncci, so we
# have to pick a different target name and override the output executable name

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(make_6j_table_q3 PRIVATE fmt::fmt am::am)

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
install(TARGETS make_6j_table_q3 spncci_exe generate_seeds_lsu3shell DESTINATION bin)
