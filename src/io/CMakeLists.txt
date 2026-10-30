add_library(arise_io
    csv_table.cpp
    number.cpp
    peak_csv.cpp)
add_library(arise::io ALIAS arise_io)

target_link_libraries(arise_io PUBLIC arise::core)
arise_set_warnings(arise_io)
