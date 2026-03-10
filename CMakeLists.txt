cmake_minimum_required(VERSION 3.30)
project(dualmc)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()

add_library(dualmc INTERFACE)
target_include_directories(dualmc INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(examples)
endif()