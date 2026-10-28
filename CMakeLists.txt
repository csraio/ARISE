cmake_minimum_required(VERSION 3.25)

project(ARISE
    VERSION 0.1.0
    DESCRIPTION "Anisotropy-Resolved Isomer Scoring Engine"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ARISE_BUILD_TESTS "Build the GoogleTest suites" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(AriseWarnings)

add_subdirectory(third-party)
add_subdirectory(src)

if(ARISE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
