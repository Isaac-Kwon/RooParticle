# 3.16 is the oldest version that works without deprecation warnings on
# cmake 3.28 (compatibility with < 3.5 is deprecated); the "..." upper bound
# opts in to the newer policy set when available.
cmake_minimum_required(VERSION 3.16...3.28)
project(trillision LANGUAGES CXX)

# by default build optimized code with debug symbols
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---------------------------------------------------------------- ROOT
# Components that the sources use (same set the old FindROOT.cmake linked).
set(RP_ROOT_COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint
                       Postscript Matrix Physics MathCore Thread MultiProc)

# Preferred: ROOTConfig.cmake (found via CMAKE_PREFIX_PATH / ROOTSYS, which
# thisroot.sh sets) giving imported targets ROOT::Core etc.
# Fallback: cmake/FindROOT.cmake driven by root-config, for ROOT installs
# without cmake config files. Force it with -DRP_ROOT_USE_MODULE=ON.
option(RP_ROOT_USE_MODULE "Use cmake/FindROOT.cmake (root-config) instead of ROOTConfig.cmake" OFF)
if(NOT RP_ROOT_USE_MODULE)
  find_package(ROOT CONFIG QUIET COMPONENTS ${RP_ROOT_COMPONENTS})
endif()
if(ROOT_FOUND)
  set(RP_ROOT_TARGETS)
  foreach(_c ${RP_ROOT_COMPONENTS})
    if(TARGET ROOT::${_c})
      list(APPEND RP_ROOT_TARGETS ROOT::${_c})
    endif()
  endforeach()
else()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  if(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
  endif()
  find_package(ROOT MODULE REQUIRED)
  # wrap the variables of the module in one imported target, so that the
  # rest of this file does not care which way ROOT was found
  add_library(RP_root_module INTERFACE)
  target_include_directories(RP_root_module SYSTEM INTERFACE ${ROOT_INCLUDE_DIRS})
  target_link_libraries(RP_root_module INTERFACE ${ROOT_LIBRARIES})
  target_compile_definitions(RP_root_module INTERFACE ${ROOT_DEFINITIONS})
  set(RP_ROOT_TARGETS RP_root_module)
endif()
message(STATUS "Using ROOT ${ROOT_VERSION} (${ROOT_DIR}${ROOT_CONFIG_EXECUTABLE})")

# ------------------------------------------------------ C++ standard
# ensure compatibility with ROOT: take the standard ROOT was built with
set(CMAKE_CXX_STANDARD 11)
if(ROOT_cxx17_FOUND OR ROOT_CXX_STANDARD EQUAL 17 OR ROOT_cxx_std_17_FOUND)
  set(CMAKE_CXX_STANDARD 17)
elseif(ROOT_cxx14_FOUND)
  set(CMAKE_CXX_STANDARD 14)
endif()
# ROOT's own imported targets may demand a higher standard (e.g. C++17/20)
# through cxx_std_NN compile features; those win over this default.
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

# check for build in the source tree
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(STATUS "Building in source tree")
endif()

# ------------------------------------------------------ git revision
# Short git SHA as the string macro RP_GIT_SHA; "unknown" without git or .git
set(RP_GIT_SHA "unknown")
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")  # .git is a file in worktrees
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE _sha
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _sha_rc
    ERROR_QUIET)
  if(_sha_rc EQUAL 0 AND NOT _sha STREQUAL "")
    set(RP_GIT_SHA "${_sha}")
  endif()
  # re-run cmake when HEAD moves (checkout); a new commit on the same branch
  # changes the ref file, which is only picked up at the next configure
  if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 "${PROJECT_SOURCE_DIR}/.git/HEAD")
  endif()
endif()
message(STATUS "Git revision: ${RP_GIT_SHA}")

# ------------------------------------------------------ library
file(GLOB RP_SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB RP_HEADERS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hpp)

add_library(rooparticle STATIC ${RP_SOURCES} ${RP_HEADERS})
target_include_directories(rooparticle PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(rooparticle PUBLIC ${RP_ROOT_TARGETS})
target_compile_definitions(rooparticle PUBLIC "RP_GIT_SHA=\"${RP_GIT_SHA}\"")
target_compile_options(rooparticle PUBLIC
  # same flags and order as before (the later -Wunused-variable wins)
  -Wall -Wpedantic -Wno-unused-variable -Wunused-variable -Wno-uninitialized)

# for now put the executables in the build root
# for compatibility with existing scripts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Data)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE rooparticle)

# ------------------------------------------------------ tests
enable_testing()
add_subdirectory(test)

# API documentation (make doc), never part of the normal build
include(cmake/Doxygen.cmake)
