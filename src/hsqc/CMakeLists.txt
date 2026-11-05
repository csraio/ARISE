add_library(arise_hsqc
    exchangeable.cpp)
add_library(arise::hsqc ALIAS arise_hsqc)

target_link_libraries(arise_hsqc PUBLIC arise::nmr1d)
arise_set_warnings(arise_hsqc)
