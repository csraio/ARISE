add_library(arise_chem
    molecule.cpp)
add_library(arise::chem ALIAS arise_chem)

target_link_libraries(arise_chem PUBLIC arise::core)
arise_set_warnings(arise_chem)
