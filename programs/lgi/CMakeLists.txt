cmake_minimum_required(VERSION 3.15)


## find Tomas' lsu3shell libraries
find_package(lsu3shell REQUIRED)

# ##############################################################################
# define executables
# ##############################################################################

add_executable(Brel_subspace_lister Brel_subspace_lister.cpp)
add_executable(get_cmf_lgi_labels get_cmf_lgi_labels.cpp)
add_executable(get_lgi_expansion get_lgi_expansion.cpp)
add_executable(get_spncci_seed_blocks get_spncci_seed_blocks.cpp)
add_executable(get_u3s_subspaces get_u3s_subspaces.cpp)
add_executable(transform_seed_blocks transform_seed_blocks.cpp)
add_executable(u3s_to_sp3rs_basis u3s_to_sp3rs_basis.cpp)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(Brel_subspace_lister fmt::fmt spncci::sp3rlib spncci::lgi)
target_link_libraries(
  get_cmf_lgi_labels
  PRIVATE fmt::fmt
          mcutils::mcutils
          spncci::lgi
          spncci::lsu3shell
          spncci::sp3rlib
          spncci::spncci
          spncci::u3shell
)

target_link_libraries(
  u3s_to_sp3rs_basis
  PRIVATE fmt::fmt
          lsu3shell::libraries
          mcutils::mcutils
          spncci::lgi
          spncci::lsu3shell
          spncci::sp3rlib
          spncci::spncci
          spncci::u3shell

)

target_link_libraries(
  get_lgi_expansion PRIVATE fmt::fmt mcutils::mcutils spncci::lgi
                            spncci::sp3rlib spncci::spncci spncci::u3shell
)
target_link_libraries(
  get_spncci_seed_blocks PRIVATE fmt::fmt mcutils::mcutils spncci::lgi
                                 spncci::sp3rlib spncci::spncci spncci::u3shell
)
target_link_libraries(get_u3s_subspaces PRIVATE am::am fmt::fmt spncci::sp3rlib)
target_link_libraries(
  transform_seed_blocks PRIVATE fmt::fmt spncci::lgi spncci::spncci
)

add_executable(su3rme_test su3rme_test.cpp)
target_link_libraries(
  su3rme_test PRIVATE 
  lsu3shell::SU3ME
  spncci::lsu3shell
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS Brel_subspace_lister get_cmf_lgi_labels get_lgi_expansion
                get_spncci_seed_blocks get_u3s_subspaces transform_seed_blocks
        DESTINATION bin
)
