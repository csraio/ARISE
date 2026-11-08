add_library(arise_filter
    exchangeable.cpp
    groups.cpp)
add_library(arise::filter ALIAS arise_filter)

# The measured suggestion windows ship inside the engine: read at configure time, so the
# build needs neither Python nor the file at run time.
set(ARISE_GROUP_SHIFTS_FILE "${PROJECT_SOURCE_DIR}/data/group_shifts.tsv")
file(READ "${ARISE_GROUP_SHIFTS_FILE}" ARISE_GROUP_SHIFTS_TEXT)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ARISE_GROUP_SHIFTS_FILE}")
configure_file(group_shifts_data.hpp.in
               "${CMAKE_CURRENT_BINARY_DIR}/generated/group_shifts_data.hpp" @ONLY)
target_include_directories(arise_filter PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")

target_link_libraries(arise_filter PUBLIC arise::chem)
arise_set_warnings(arise_filter)
