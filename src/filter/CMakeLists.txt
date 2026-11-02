add_library(arise_filter
    groups.cpp)
add_library(arise::filter ALIAS arise_filter)

target_link_libraries(arise_filter PUBLIC arise::chem)
arise_set_warnings(arise_filter)
