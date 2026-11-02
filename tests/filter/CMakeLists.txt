add_executable(arise_filter_tests
    groups_test.cpp
    holdout_test.cpp)

target_link_libraries(arise_filter_tests
    PRIVATE
        arise::filter
        arise_group_shifts
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_filter_tests)

gtest_discover_tests(arise_filter_tests)
