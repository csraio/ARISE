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

# Golden rule 4: never fuse a*b+c into one FMA. Clang does so by default on arm64 when
# optimizing, but not on x86-64 or wasm, so results would differ in the last bit.
# Guarded by tests/core/platform_test.cpp.
if(NOT MSVC)
    add_compile_options(-ffp-contract=off)
endif()

add_subdirectory(third-party)
add_subdirectory(src)

if(ARISE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
