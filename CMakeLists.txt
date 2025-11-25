cmake_minimum_required(VERSION 3.28)
project(dualmc)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()

add_library(dualmc INTERFACE)
target_include_directories(dualmc INTERFACE "${CMAKE_SOURCE_DIR}/include")

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(examples)
endif()