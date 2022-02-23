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
)

add_executable(generate_lgi_expansion_lsu3shell generate_lgi_expansion_lsu3shell.cpp)
target_link_libraries(
  generate_lgi_expansion_lsu3shell PRIVATE
  fmt::fmt
  am::am
  spncci::utilities
  spncci::lgi
)

# add_executable(generate_seeds_lsu3shell generate_seeds_lsu3shell.cpp)
# target_link_libraries(
#   generate_seeds_lsu3shell PRIVATE
#   fmt::fmt
#   am::am
#   spncci::utilities
#   spncci::lgi
# )


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



# add_executable(Brel_subspace_lister EXCLUDE_FROM_ALL Brel_subspace_lister.cpp)
# target_link_libraries(Brel_subspace_lister fmt::fmt spncci::sp3rlib spncci::lgi)
# add_executable(get_cmf_lgi_labels EXCLUDE_FROM_ALL get_cmf_lgi_labels.cpp)
# target_link_libraries(
#   get_cmf_lgi_labels
#   PRIVATE fmt::fmt
#           mcutils::mcutils
#           spncci::lgi
#           spncci::lsu3shell
#           spncci::sp3rlib
#           spncci::spncci
#           spncci::u3shell
# )
# add_executable(u3s_to_sp3rs_basis EXCLUDE_FROM_ALL u3s_to_sp3rs_basis.cpp)
# target_link_libraries(
#   u3s_to_sp3rs_basis
#   PRIVATE fmt::fmt
#           lsu3shell::libraries
#           mcutils::mcutils
#           spncci::lgi
#           spncci::lsu3shell
#           spncci::sp3rlib
#           spncci::spncci
#           spncci::u3shell
# )

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS generate_lgi_expansion_lsu3shell
        DESTINATION bin
)
