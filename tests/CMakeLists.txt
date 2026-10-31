include(GoogleTest)

# Recorded real-world inputs, embedded at configure time so tests need no network or
# file access (also under wasm). Link `arise_test_fixtures` to use them.
set(ARISE_PUBCHEM_FIXTURE "${CMAKE_CURRENT_SOURCE_DIR}/data/pubchem_C4H8O2.csv")
file(READ "${ARISE_PUBCHEM_FIXTURE}" ARISE_PUBCHEM_C4H8O2)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ARISE_PUBCHEM_FIXTURE}")
configure_file(fixtures/pubchem_fixture.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/fixtures/pubchem_fixture.hpp" @ONLY)
add_library(arise_test_fixtures INTERFACE)
target_include_directories(arise_test_fixtures INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/fixtures")

add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(nmr1d)
add_subdirectory(chem)
