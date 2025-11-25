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

add_example(atom_example)
add_example(atom_example_visual)

#[[
set(EXAMPLE_APP_SOURCES
    apps/example/example.cpp
    apps/example/main.cpp
)

set(GENTABLES_APP_SOURCES
    apps/gentables/gentables.cpp
    apps/gentables/main.cpp
)

# build application
add_executable(dmc ${EXAMPLE_APP_SOURCES})
add_executable(gentables ${GENTABLES_APP_SOURCES})
]]