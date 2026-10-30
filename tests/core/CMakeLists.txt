add_executable(arise_core_tests
    units_test.cpp
    result_test.cpp
    platform_test.cpp)

target_link_libraries(arise_core_tests
    PRIVATE
        arise::core
        GTest::gtest_main)
arise_set_warnings(arise_core_tests)

gtest_discover_tests(arise_core_tests)
