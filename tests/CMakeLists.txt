include(GoogleTest)

# Recorded real-world inputs, embedded at configure time so tests need no network or
# file access (also under wasm). Link `arise_test_fixtures` to use them.
set(ARISE_PUBCHEM_FIXTURE "${CMAKE_CURRENT_SOURCE_DIR}/data/pubchem_C4H8O2.csv")
file(READ "${ARISE_PUBCHEM_FIXTURE}" ARISE_PUBCHEM_C4H8O2)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ARISE_PUBCHEM_FIXTURE}")
configure_file(fixtures/pubchem_fixture.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/fixtures/pubchem_fixture.hpp" @ONLY)
# Reference data is larger than the 65,536 characters a string literal may portably hold,
# so it is embedded as raw-string chunks that tests join at run time. The ""sv suffix gives
# each chunk's length directly; constructing a string_view from a plain literal would count
# 60,000 characters at compile time and exceed clang's constexpr step limit with libstdc++.
function(arise_embed_chunks file out_chunks out_count)
    file(READ "${file}" content)
    string(LENGTH "${content}" length)
    set(chunk_size 60000)
    set(chunks "")
    set(count 0)
    set(offset 0)
    while(offset LESS length)
        string(SUBSTRING "${content}" ${offset} ${chunk_size} piece)
        string(APPEND chunks "    R\"ref(${piece})ref\"sv,\n")
        math(EXPR offset "${offset} + ${chunk_size}")
        math(EXPR count "${count} + 1")
    endwhile()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${file}")
    set(${out_chunks} "${chunks}" PARENT_SCOPE)
    set(${out_count} ${count} PARENT_SCOPE)
endfunction()

set(ARISE_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/reference")
arise_embed_chunks("${ARISE_REFERENCE_DIR}/nmrshiftdb2_13c_sample.tsv"
                   ARISE_NMRSHIFTDB_CHUNKS ARISE_NMRSHIFTDB_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/pubchem_rdkit_oracle.tsv"
                   ARISE_ORACLE_CHUNKS ARISE_ORACLE_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/chebi_classes.tsv"
                   ARISE_CHEBI_CHUNKS ARISE_CHEBI_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/known_differences.tsv"
                   ARISE_KNOWN_DIFFERENCES_CHUNKS ARISE_KNOWN_DIFFERENCES_CHUNK_COUNT)
configure_file(fixtures/reference_fixtures.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/fixtures/reference_fixtures.hpp" @ONLY)
add_library(arise_test_fixtures INTERFACE)
target_include_directories(arise_test_fixtures INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/fixtures")

add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(nmr1d)
add_subdirectory(chem)
