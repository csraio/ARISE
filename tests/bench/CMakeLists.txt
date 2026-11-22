add_executable(arise_bench_tests prune_pipeline_test.cpp snapshot_test.cpp)
target_link_libraries(arise_bench_tests
    PRIVATE
        arise_bench_pipeline
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_bench_tests)

gtest_discover_tests(arise_bench_tests)
