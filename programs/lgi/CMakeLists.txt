cmake_minimum_required(VERSION 3.15)

###############################################################################
#define executables and link dependencies
###############################################################################
if(BUILD_LSU3SHELL)
  message(STATUS "Including program/lgi/generate_lgi_expansion_lsu3shell in the build.")
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
  message(STATUS "Including program/lgi/get_u3s_subspaces in the build.")
  add_executable(get_u3s_subspaces EXCLUDE_FROM_ALL get_u3s_subspaces.cpp)
  target_link_libraries(get_u3s_subspaces PRIVATE am::am fmt::fmt spncci::sp3rlib spncci::u3ncsm spncci::utilities)
endif()


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
if(BUILD_LSU3SHELL)
  install(TARGETS generate_lgi_expansion_lsu3shell  get_u3s_subspaces DESTINATION bin)
endif()
