cmake_minimum_required(VERSION 3.20)

# add headers and sources
set(su4lib_headers u4_branching.h 
)
# set(su4lib_sources )

# define am library and add source files
add_library(su4lib INTERFACE)
add_library(su4lib::su4lib ALIAS su4lib)

# ##############################################################################
# link dependencies
# ##############################################################################



# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH su4lib_PARENT_DIR)
target_include_directories(
  su4lib INTERFACE $<BUILD_INTERFACE:${su4lib_PARENT_DIR}>
                    $<INSTALL_INTERFACE:include/spncci>
)

target_link_libraries(
  su4lib INTERFACE am::am
         # fmt::fmt
)
# ##############################################################################
# define installed targets
#
# actual installation is performed by top-level CMakeLists.txt
# ##############################################################################
install(
  TARGETS su4lib
  DESTINATION lib
  EXPORT spncciTargets
)
install(FILES ${su4lib_headers} DESTINATION include/spncci/su4lib)

# ##############################################################################
# define tests
# ##############################################################################

set(su4lib_TEST
    u4_branching_test
)

add_custom_target(su4lib_tests)
foreach(test_name IN LISTS su4lib_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
  target_link_libraries(${test_name} su4lib::su4lib)
  add_dependencies(su4lib_tests ${test_name})
endforeach()

target_link_libraries(u4_branching_test fmt::fmt Eigen3::Eigen)


if(SPNCCI_MASTER_PROJECT)
  add_dependencies(tests su4lib_tests)
endif()
