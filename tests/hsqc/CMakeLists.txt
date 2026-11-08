add_executable(arise_hsqc_tests
    bmrb_test.cpp
    exchangeable_test.cpp
    nmrshiftdb2_test.cpp)

# arise::chem is linked for the BMRB check only, to tell which BMRB atoms are carbons; the
# hsqc module itself does not depend on it.
target_link_libraries(arise_hsqc_tests
    PRIVATE
        arise::hsqc
        arise::chem
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_hsqc_tests)

gtest_discover_tests(arise_hsqc_tests)
