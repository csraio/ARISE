add_executable(arise_nmr1d_tests
    dept_test.cpp
    solvent_test.cpp)

target_link_libraries(arise_nmr1d_tests
    PRIVATE
        arise::nmr1d
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_nmr1d_tests)

gtest_discover_tests(arise_nmr1d_tests)
