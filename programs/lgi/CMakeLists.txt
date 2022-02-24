cmake_minimum_required(VERSION 3.15)


## find Tomas' lsu3shell libraries
find_package(lsu3shell REQUIRED)

###############################################################################
#define executables and link dependencies
###############################################################################

add_executable(su3rme_test su3rme_test.cpp)
target_link_libraries(
  su3rme_test PRIVATE
  lsu3shell::SU3ME
  spncci::lsu3shell
  spncci::u3ncsm
)

add_executable(generate_lgi_expansion_lsu3shell generate_lgi_expansion_lsu3shell.cpp)
target_link_libraries(
  generate_lgi_expansion_lsu3shell PRIVATE
  fmt::fmt
  am::am
  spncci::utilities
  spncci::spncci_basis
  spncci::u3ncsm
  spncci::lgi
)


## For lanczos decompositions
add_executable(get_u3s_subspaces EXCLUDE_FROM_ALL get_u3s_subspaces.cpp)
target_link_libraries(get_u3s_subspaces PRIVATE am::am fmt::fmt spncci::sp3rlib)

## Old version of generating seeds
add_executable(get_lgi_expansion EXCLUDE_FROM_ALL get_lgi_expansion.cpp)
target_link_libraries(get_lgi_expansion PRIVATE fmt::fmt mcutils::mcutils spncci::lgi
                      spncci::sp3rlib spncci::spncci spncci::u3shell)

add_executable(get_spncci_seed_blocks EXCLUDE_FROM_ALL get_spncci_seed_blocks.cpp)
target_link_libraries(get_spncci_seed_blocks PRIVATE fmt::fmt mcutils::mcutils spncci::lgi
                        spncci::sp3rlib spncci::spncci spncci::u3shell)


add_executable(transform_seed_blocks EXCLUDE_FROM_ALL transform_seed_blocks.cpp)
target_link_libraries(transform_seed_blocks PRIVATE fmt::fmt spncci::lgi spncci::spncci)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS generate_lgi_expansion_lsu3shell
        DESTINATION bin
)
