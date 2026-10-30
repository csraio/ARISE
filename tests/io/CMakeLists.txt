add_executable(arise_io_tests
    carbon_csv_test.cpp
    number_test.cpp
    peak_csv_test.cpp)

target_link_libraries(arise_io_tests
    PRIVATE
        arise::io
        GTest::gtest_main)
arise_set_warnings(arise_io_tests)

gtest_discover_tests(arise_io_tests)
