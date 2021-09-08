cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# define executables
# ##############################################################################

add_executable(
  generate_relative_u3st_operators generate_relative_u3st_operators.cpp
)
add_executable(truncate_interaction truncate_interaction.cpp)

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(
  generate_relative_u3st_operators PRIVATE Boost::headers fmt::fmt
                                           mcutils::mcutils spncci::u3shell
)
target_link_libraries(truncate_interaction PRIVATE fmt::fmt spncci::u3shell)

# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(TARGETS generate_relative_u3st_operators truncate_interaction
        DESTINATION bin
)
