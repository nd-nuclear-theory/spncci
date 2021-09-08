cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

add_executable(spectra_dense_timing spectra_dense_timing.cpp)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  spectra_dense_timing Eigen3::Eigen mcutils::mcutils Spectra::Spectra
)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS spectra_dense_timing DESTINATION bin)
