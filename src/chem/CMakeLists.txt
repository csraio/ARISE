add_library(arise_chem
    molecule.cpp
    perception.cpp)
add_library(arise::chem ALIAS arise_chem)

target_link_libraries(arise_chem PUBLIC arise::core arise::nmr1d)
arise_set_warnings(arise_chem)
