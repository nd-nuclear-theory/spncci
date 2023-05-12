###############################################################################
# FetchContent declarations for external libraries
#
# Anna E. McCoy
# Institute for Nuclear Theory 
#
# SPDX-License-Identifier: MIT
#
# 3/29/22 (aem) : Created
###############################################################################
cmake_minimum_required(VERSION 3.20)

include(FetchContent)
FetchContent_Declare(
  wigxjpf
  GIT_REPOSITORY https://github.com/nd-nuclear-theory/wigxjpf.git
  GIT_TAG        main
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  ompilancz
  GIT_REPOSITORY https://gitlab.com/nd-nuclear-theory/ompilancz.git
  GIT_TAG        c7ee7d27e1a38ad2581600a63b8265e3ad3fde22
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  SU3lib
  GIT_REPOSITORY https://gitlab.com/nd-nuclear-theory/SU3lib.git
  GIT_TAG        master
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  su3wrc
  GIT_REPOSITORY https://github.com/nd-nuclear-theory/su3wrc.git
  GIT_TAG        cmake
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  ndsu3lib
  GIT_REPOSITORY https://github.com/nd-nuclear-theory/ndsu3lib.git
  GIT_TAG        master
  GIT_SHALLOW    TRUE
)


FetchContent_Declare(
  GTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  lsu3shell
#  GIT_REPOSITORY https://github.com/nd-nuclear-theory/lsu3shell.git
  GIT_REPOSITORY https://gitlab.com/nd-nuclear-theory/lsu3shell.git
#  GIT_REPOSITORY git@gitlab.com:nd-nuclear-theory/lsu3shell.git
  GIT_TAG        master
  GIT_SHALLOW    TRUE
)


