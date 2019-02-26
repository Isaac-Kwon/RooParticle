cmake_minimum_required(VERSION 3.2)
project(trillision)

# by default build optimized code with debug symbols
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# look for ROOT and warn if unavailable
message(STATUS "Looking for ROOT")
# try first with ROOTConfig.cmake,
# then fallback to FindROOT.cmake (non-cmake installation of ROOT)
if(DISABLE_ROOT_CONFIG STREQUAL "OFF")
  find_package(ROOT CONFIG QUIET)
endif()
if(NOT ROOT_FOUND)
  list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(ROOT MODULE QUIET REQUIRED)
endif()
if(ROOT_FOUND)
  message(STATUS "Using ROOT ${ROOT_VERSION}")
else()
  message(WARNING "Compilation without ROOT results in limited functionality!")
endif()

# require C++11-capable compiler
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
    message(WARNING "Compiler with incomplete C++11 support detected: functionality limited")
    set(OLD_CXX true)
  endif()
endif()

# set C++ standard ensuring compatibility with ROOT
set(CMAKE_CXX_STANDARD 11)
if(ROOT_cxx17_FOUND)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
elseif(ROOT_cxx14_FOUND)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
elseif(ROOT_cxx11_FOUND)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

# warn about undefined symbols already for libraries
#link_libraries("-Wl,--as-needed")
link_libraries("-Wl,${ROOT_LIBRARIES}")
# link_libraries("${ROOT_EXE_LINKER_FLAGS},${ROOT_LIBRARIES}")
# link_libraries("-Wl")
add_compile_options("-Wall")
add_compile_options("-Werror")
add_compile_options("-Wpedantic")
#add_compile_options("-v")

# check for build in the source tree
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(STATUS "Building in source tree")
  set(BUILDINSOURCE true)
endif()

message(STATUS "In Project Directory: ${PROJECT_SOURCE_DIR}")
message(STATUS "ROOT Include Directory: ${ROOT_INCLUDE_DIRS}")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)

# add subdirectories for the actual work
#include_directories(src)
include_directories(include)
include_directories(${ROOT_INCLUDE_DIRS})

add_subdirectory(test)

add_executable(main main.cpp ${sources} ${headers})

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Data)
