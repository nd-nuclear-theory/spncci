cmake_minimum_required(VERSION 3.15)

# ##############################################################################
# internal dependencies
#
# note: since this script is include-d from the parent directory, we must
# prepend "programs/" to each path
# ##############################################################################

add_subdirectory(programs/lgi)
add_subdirectory(programs/linear_algebra)
add_subdirectory(programs/operators)
add_subdirectory(programs/rotor)
add_subdirectory(programs/spncci)
add_subdirectory(programs/su3calc)
add_subdirectory(programs/unit_tensors)
