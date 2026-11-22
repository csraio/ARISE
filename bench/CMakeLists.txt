# Development-only benchmark (plan 3d): never part of the engine, CLI or web app. The
# benchmark data are embedded through the test fixtures, so it builds with the tests.
add_library(arise_bench_pipeline STATIC prune_pipeline.cpp snapshot.cpp)
target_include_directories(arise_bench_pipeline PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(arise_bench_pipeline PUBLIC arise::filter arise::io)
arise_set_warnings(arise_bench_pipeline)

add_executable(arise_bench arise_bench.cpp)
target_link_libraries(arise_bench PRIVATE arise_bench_pipeline arise_test_fixtures)
arise_set_warnings(arise_bench)

if(NOT EMSCRIPTEN)
    add_executable(arise_snapshot arise_snapshot.cpp)
    target_link_libraries(arise_snapshot PRIVATE arise_bench_pipeline arise_test_fixtures)
    arise_set_warnings(arise_snapshot)
endif()
