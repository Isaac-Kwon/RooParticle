# Every test/main_<name>.cpp becomes the executable test_<name>.
file(GLOB RP_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/main_*.cpp)

# Files known not to build. Delete a line as soon as the file compiles.
set(RP_TEST_EXCLUDE
  # (empty: add "main_<name>.cpp  # reason" here)
  # does not compile: undeclared randomseed/nparticle/imp_min/imp_max, stale inspector ctor, TString->Float_t
  main_rutherford_track.cpp
)

set(RP_TEST_TARGETS)
foreach(src ${RP_TEST_SOURCES})
  get_filename_component(file ${src} NAME)
  if(file IN_LIST RP_TEST_EXCLUDE)
    message(STATUS "Test excluded: ${file}")
    continue()
  endif()
  get_filename_component(name ${file} NAME_WE)
  string(SUBSTRING ${name} 5 -1 name)  # strip "main_"
  set(exe test_${name})
  add_executable(${exe} ${src})
  target_link_libraries(${exe} PRIVATE rooparticle)
  list(APPEND RP_TEST_TARGETS ${exe})

  # Only self-checking tests (no arguments, non-zero exit on failure) go to
  # ctest. The other executables are long simulations or open GUI canvases.
  # They write to Data/ relative to the working directory, hence the
  # WORKING_DIRECTORY.
  if(name STREQUAL "integrator_compare" OR name MATCHES "_check$")
    add_test(NAME ${exe} COMMAND ${exe} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endif()
endforeach()

install(TARGETS ${RP_TEST_TARGETS} RUNTIME DESTINATION bin)
