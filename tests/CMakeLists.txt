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
arise_embed_chunks("${ARISE_REFERENCE_DIR}/measurements.tsv"
                   ARISE_MEASUREMENTS_CHUNKS ARISE_MEASUREMENTS_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/spec_rdkit_oracle.tsv"
                   ARISE_SPEC_ORACLE_CHUNKS ARISE_SPEC_ORACLE_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/valence_probes.tsv"
                   ARISE_VALENCE_PROBES_CHUNKS ARISE_VALENCE_PROBES_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/dou_checks.tsv"
                   ARISE_DOU_CHECKS_CHUNKS ARISE_DOU_CHECKS_CHUNK_COUNT)
# The signed-off definitions, so tests can hold the doc's ChEBI IDs and numbers to the data.
arise_embed_chunks("${PROJECT_SOURCE_DIR}/docs/chem_group_definitions.md"
                   ARISE_GROUP_DEFINITIONS_CHUNKS ARISE_GROUP_DEFINITIONS_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/nmrshiftdb2_13c_all.tsv"
                   ARISE_NMRSHIFTDB_ALL_CHUNKS ARISE_NMRSHIFTDB_ALL_CHUNK_COUNT)
arise_embed_chunks("${PROJECT_SOURCE_DIR}/data/group_shifts.tsv"
                   ARISE_GROUP_SHIFTS_CHUNKS ARISE_GROUP_SHIFTS_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/pubchem_cleaning_oracle.tsv"
                   ARISE_CLEANING_ORACLE_CHUNKS ARISE_CLEANING_ORACLE_CHUNK_COUNT)
arise_embed_chunks("${ARISE_REFERENCE_DIR}/named_compounds.tsv"
                   ARISE_NAMED_COMPOUNDS_CHUNKS ARISE_NAMED_COMPOUNDS_CHUNK_COUNT)
# Every recorded PubChem formula search, each file's chunks preceded by a "## <formula>"
# chunk, for the cleaning check.
set(ARISE_PUBCHEM_SEARCHES_CHUNKS "")
set(ARISE_PUBCHEM_SEARCHES_CHUNK_COUNT 0)
file(GLOB ARISE_PUBCHEM_SEARCHES "${CMAKE_CURRENT_SOURCE_DIR}/data/pubchem_*.csv")
list(SORT ARISE_PUBCHEM_SEARCHES)
foreach(search ${ARISE_PUBCHEM_SEARCHES})
    get_filename_component(stem "${search}" NAME_WE)
    string(REPLACE "pubchem_" "" formula "${stem}")
    arise_embed_chunks("${search}" chunks count)
    string(APPEND ARISE_PUBCHEM_SEARCHES_CHUNKS "    R\"ref(## ${formula}\n)ref\"sv,\n${chunks}")
    math(EXPR ARISE_PUBCHEM_SEARCHES_CHUNK_COUNT
         "${ARISE_PUBCHEM_SEARCHES_CHUNK_COUNT} + ${count} + 1")
endforeach()

# The hand-written tests, so a test can require RDKit's view of every SMILES in them.
set(ARISE_SPEC_SOURCES_CHUNKS "")
set(ARISE_SPEC_SOURCES_CHUNK_COUNT 0)
foreach(source chem/molecule_test.cpp chem/perception_test.cpp filter/groups_test.cpp)
    arise_embed_chunks("${CMAKE_CURRENT_SOURCE_DIR}/${source}" chunks count)
    string(APPEND ARISE_SPEC_SOURCES_CHUNKS "${chunks}")
    math(EXPR ARISE_SPEC_SOURCES_CHUNK_COUNT "${ARISE_SPEC_SOURCES_CHUNK_COUNT} + ${count}")
endforeach()
configure_file(fixtures/reference_fixtures.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/fixtures/reference_fixtures.hpp" @ONLY)
add_library(arise_test_fixtures INTERFACE)
target_include_directories(arise_test_fixtures INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/fixtures")

add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(nmr1d)
add_subdirectory(chem)
add_subdirectory(filter)
