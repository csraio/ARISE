add_executable(arise_io_tests
    candidates_test.cpp
    carbon_csv_test.cpp
    correlation_csv_test.cpp
    number_test.cpp
    peak_csv_test.cpp
    pubchem_test.cpp)

target_link_libraries(arise_io_tests
    PRIVATE
        arise::io
        GTest::gtest_main)
arise_set_warnings(arise_io_tests)

# Embed recorded responses so tests need no network or file access (also under wasm).
set(ARISE_PUBCHEM_FIXTURE "${CMAKE_CURRENT_SOURCE_DIR}/data/pubchem_C4H8O2.csv")
file(READ "${ARISE_PUBCHEM_FIXTURE}" ARISE_PUBCHEM_C4H8O2)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ARISE_PUBCHEM_FIXTURE}")
configure_file(pubchem_fixture.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/pubchem_fixture.hpp" @ONLY)
target_include_directories(arise_io_tests PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

gtest_discover_tests(arise_io_tests)
