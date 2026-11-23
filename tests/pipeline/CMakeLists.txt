add_executable(arise_pipeline_tests
    json_test.cpp
    pipeline_test.cpp)

target_link_libraries(arise_pipeline_tests
    PRIVATE
        arise::pipeline
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_pipeline_tests)

gtest_discover_tests(arise_pipeline_tests)
