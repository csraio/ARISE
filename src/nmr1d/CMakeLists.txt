add_library(arise_nmr1d
    dept.cpp
    solvent.cpp)
add_library(arise::nmr1d ALIAS arise_nmr1d)

# The residual-solvent and water peaks ship inside the engine: read at configure time, so
# the build needs neither Python nor the file at run time.
set(ARISE_SOLVENT_PEAKS_FILE "${PROJECT_SOURCE_DIR}/data/solvent_peaks.tsv")
file(READ "${ARISE_SOLVENT_PEAKS_FILE}" ARISE_SOLVENT_PEAKS_TEXT)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ARISE_SOLVENT_PEAKS_FILE}")
configure_file(solvent_peaks_data.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/generated/solvent_peaks_data.hpp" @ONLY)
target_include_directories(arise_nmr1d PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")

target_link_libraries(arise_nmr1d PUBLIC arise::io)
arise_set_warnings(arise_nmr1d)
