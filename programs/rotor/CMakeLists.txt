cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

# add_executable(sp3r_rotor sp3r_rotor.cpp)
add_executable(transitions transitions.cpp)
add_executable(triaxial_rotor triaxial_rotor.cpp)

# ##############################################################################
# link dependencies
# ##############################################################################

# target_link_libraries(
#   sp3r_rotor
#   PRIVATE Boost::headers
#           Eigen3::Eigen
#           fmt::fmt
#           am::am
#           m
#           spncci::lgi
#           spncci::sp3rlib
# )
target_link_libraries(
  transitions
  PRIVATE Boost::headers
          m
          fmt::fmt
          am::am
          mcutils::mcutils
          spncci::lgi
          spncci::sp3rlib
          Eigen3::Eigen
)
target_link_libraries(
  triaxial_rotor PRIVATE m Eigen3::Eigen am::am spncci::sp3rlib
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS
  triaxial_rotor
  # sp3r_rotor
  transitions
  DESTINATION bin
  )
