# Development-only tools (built with the tests, never part of the engine, CLI or web app).
add_library(arise_group_shifts STATIC group_shifts.cpp)
target_include_directories(arise_group_shifts PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(arise_group_shifts PUBLIC arise::chem)
arise_set_warnings(arise_group_shifts)

if(NOT EMSCRIPTEN)
    add_executable(arise_measure_group_shifts measure_group_shifts.cpp)
    target_link_libraries(arise_measure_group_shifts PRIVATE arise_group_shifts)
    arise_set_warnings(arise_measure_group_shifts)
endif()
