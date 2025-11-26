function(ADD_EXAMPLE name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT ARG_SOURCES)
        set(ARG_SOURCES "${name}.cpp")
    endif()

    add_executable(${name} ${ARG_SOURCES})
    target_link_libraries(${name} PRIVATE dualmc)
    if(ARG_LIBS)
        target_link_libraries(${name} PRIVATE ${ARG_LIBS})
    endif()

    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_LIST_DIR})

    set_target_properties(${name} PROPERTIES FOLDER "Examples")
    if (MSVC)
        target_compile_options(${name} PRIVATE /W4)
    else()
        target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic -Wno-missing-field-initializers)
    endif()
endfunction()

file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.0/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

add_example(atom_example)
add_example(atom_example_visual)