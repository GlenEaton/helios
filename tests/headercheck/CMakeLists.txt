if(HELIOS_HEADERCHECK)
  file(GLOB_RECURSE helios_headers
    "${PROJECT_SOURCE_DIR}/src/*.h"
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
  )

  foreach(header ${helios_headers})
    # Check that this is not from the pybinds subdirectory. We excllude that directory
    # as it would complicate this logic to link properly against the Python standard library.
    string(FIND "${header}" "pybinds" pybind_index)
    string(FIND "${header}" "src/test" test_index)
    if((pybind_index EQUAL -1) AND (test_index EQUAL -1))
      get_filename_component(header_name ${header} NAME)
      string(REPLACE ".h" "" header_name ${header_name})
      string(REPLACE ".hpp" "" header_name ${header_name})
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/template.cpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/${header_name}.headercheck.cpp
        @ONLY
      )
      add_executable(${header_name}_headercheck
        ${CMAKE_CURRENT_BINARY_DIR}/${header_name}.headercheck.cpp
      )
      target_link_libraries(${header_name}_headercheck PUBLIC helios)
    endif()
  endforeach()
endif()
