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
    target_compile_features(${name} PRIVATE cxx_std_23)

    set_target_properties(${name} PROPERTIES FOLDER "Examples")
    if (MSVC)
        target_compile_options(${name} PRIVATE /W4)
    else()
        target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic -Wno-missing-field-initializers)
    endif()

    if(WIN32)
        add_custom_command(TARGET ${name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:${name}> $<TARGET_FILE_DIR:${name}>
            COMMAND_EXPAND_LISTS
        )
    endif()
endfunction()

file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.0/CPM.cmake
  ${CMAKE_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a
)
include(${CMAKE_BINARY_DIR}/cmake/CPM.cmake)

CPMAddPackage(
  NAME daxa
  GITHUB_REPOSITORY ShatteredWorld/Daxa
  GIT_TAG bb492aa7e2d68e4e784b3e79bb1f501cb22876e7
  OPTIONS
    "DAXA_ENABLE_UTILS_FSR3 OFF"
    "DAXA_ENABLE_UTILS_IMGUI OFF"
    "DAXA_ENABLE_UTILS_MEM ON"
    "DAXA_ENABLE_UTILS_PIPELINE_MANAGER_GLSLANG OFF"
    "DAXA_ENABLE_UTILS_PIPELINE_MANAGER_SLANG ON"
    "DAXA_ENABLE_UTILS_PIPELINE_MANAGER_SPIRV_VALIDATION OFF"
    "DAXA_ENABLE_UTILS_TASK_GRAPH ON"
    "DAXA_ENABLE_TESTS OFF"
    "DAXA_ENABLE_TOOLS OFF"
    "DAXA_ENABLE_STATIC_ANALYSIS OFF"
)

CPMAddPackage(
  NAME glm
  GITHUB_REPOSITORY g-truc/glm
  GIT_TAG 1.0.1
)

CPMAddPackage(
  NAME glfw
  GITHUB_REPOSITORY glfw/glfw
  GIT_TAG 3.4
  OPTIONS
    "GLFW_BUILD_TESTS OFF"
    "GLFW_BUILD_EXAMPLES OFF"
    "GLFW_BUILD_DOCS OFF"
)

CPMAddPackage(
    NAME FastNoiseLite
    GITHUB_REPOSITORY Auburn/FastNoiseLite
    GIT_TAG v1.1.1
    DOWNLOAD_ONLY YES
)
if(FastNoiseLite_ADDED)
  add_library(FastNoiseLite INTERFACE IMPORTED)
  target_include_directories(FastNoiseLite INTERFACE "${FastNoiseLite_SOURCE_DIR}/Cpp")
endif()

add_example(atom_example)
add_example(atom_example_visual LIBS daxa glm glfw FastNoiseLite)
