
set(TEST_EXE_ROOT
  test_print
  test_rutherford
  test_collider_HeHe
  test_3body_uniform
  test_3body_approx
)

# for now put the executables in the build root
# for compatibility with existing scripts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

foreach(exe ${TEST_EXE_ROOT})
  string(REGEX REPLACE "^test_" "" test ${exe})
  add_executable(${exe} main_${test}.cpp ${sources} ${headers})
  target_include_directories(${exe} PRIVATE ${ROOT_INCLUDE_DIRS})
  #target_link_libraries(${exe} alpide)
  #include_directories(include)

endforeach()

#include_directories(include)

install(TARGETS
  ${TEST_EXE_ROOT}
  RUNTIME DESTINATION bin
)
