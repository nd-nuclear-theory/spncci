Notes for installing spncci with cmake at NERSC

If linking to lsu3shell, required libraries are:
  SU3lib: https://gitlab.com/nd-nuclear-theory/su3lib
  ompilancz:  https://gitlab.com/nd-nuclear-theory/ompilancz
  lsu3shell: https://gitlab.com/nd-nuclear-theory/lsu3shell

After compiling libraries, set environment variables for location of 
	SU3lib_DIR: path to SU3libConfig.cmake
	ompilancz_DIR: path to ompilanczConfig.cmake
  lsu3shell_DIR: path to lsu3shell build directory

If environment variables not set, then cmake will try to clone and build as part of spncci build.

For compiling on NERSC: set environment variables in env-intel-nersc.sh and source to load requisite modules.

env-intel-nersc.sh will also set some require environment variables

Compiling:
cmake -B <build dir>

	To build with debug turned on	
	cmake -DCMAKE_BUILD_TYPE=Debug

	To build with debug off.  Note, sets -DNDEBUG which turns off asserts
	cmake -DCMAKE_BUILD_TYPE=Release

cmake --build <build dir> -j <N>

Example: Set up build directory with linking to LSU3Shell turned off:
  cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug -DBUILD_LSU3SHELL=OFF

Example: Build specific target
  cmake --build <build dir> -j <N> --target recurrence_indexing_test


To change compile options, e.g., build with lsu3shell option, run ccmake on build
directory.

------------------------------------------------------------------------------------
cmake 3.22 bug: Not correctly identifying Cray Programming Environment.  When
setting up build, explicitly set compilers:

  env CC=cc CXX=CC FC=ftn cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug

See https://docs.nersc.gov/development/build-tools/cmake/#use-the-cray-compiler-wrappers.

Should see near begining of the output:
  - Cray Programming Environment <version and compiler>
