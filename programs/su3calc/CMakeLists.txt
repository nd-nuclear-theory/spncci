cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

add_executable(sp3rcalc sp3rcalc.cpp)
add_executable(su3calc su3calc.cpp)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  sp3rcalc PRIVATE fmt::fmt mcutils::mcutils spncci::sp3rlib
)
target_link_libraries(su3calc PRIVATE Eigen3::Eigen fmt::fmt spncci::sp3rlib)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS sp3rcalc su3calc DESTINATION bin)
