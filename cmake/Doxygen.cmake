# Optional API documentation. Add `include(cmake/Doxygen.cmake)` to CMakeLists.txt.
# Creates the target `doc` (not part of `all`); never makes the normal build fail.

get_filename_component(RP_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

file(GLOB _rp_doxygen_hints "$ENV{HOME}/opt/doxygen-*/bin")
find_program(DOXYGEN_EXECUTABLE doxygen HINTS ${_rp_doxygen_hints})
find_package(Doxygen QUIET COMPONENTS doxygen OPTIONAL_COMPONENTS dot)

if(DOXYGEN_FOUND)
  set(RP_DOC_OUTPUT_DIR "${CMAKE_BINARY_DIR}/docs")
  if(DOXYGEN_DOT_FOUND)
    set(RP_HAVE_DOT YES)
  else()
    set(RP_HAVE_DOT NO)
  endif()
  configure_file("${RP_SOURCE_DIR}/docs/Doxyfile.in" "${RP_DOC_OUTPUT_DIR}/Doxyfile" @ONLY)
  # no dependency on other targets, so it works wherever it is included
  add_custom_target(doc
    COMMAND ${CMAKE_COMMAND} -E make_directory "${RP_DOC_OUTPUT_DIR}"
    COMMAND ${DOXYGEN_EXECUTABLE} "${RP_DOC_OUTPUT_DIR}/Doxyfile"
    COMMAND ${CMAKE_COMMAND} -E echo "Docs: ${RP_DOC_OUTPUT_DIR}/html/index.html (warnings: ${RP_DOC_OUTPUT_DIR}/doxygen-warnings.log)"
    WORKING_DIRECTORY "${RP_SOURCE_DIR}"
    COMMENT "Generating API documentation with Doxygen ${DOXYGEN_VERSION}"
    VERBATIM)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${RP_SOURCE_DIR}/docs/Doxyfile.in")
  message(STATUS "Doxygen ${DOXYGEN_VERSION} found: target 'doc' available (graphviz dot: ${RP_HAVE_DOT})")
else()
  message(STATUS "Doxygen not found: target 'doc' not available (normal build unaffected)")
endif()
