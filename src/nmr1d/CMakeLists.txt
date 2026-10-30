add_library(arise_nmr1d
    dept.cpp)
add_library(arise::nmr1d ALIAS arise_nmr1d)

target_link_libraries(arise_nmr1d PUBLIC arise::io)
arise_set_warnings(arise_nmr1d)
