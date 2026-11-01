add_executable(arise_chem_tests
    molecule_test.cpp
    perception_test.cpp
    reference_test.cpp)

target_link_libraries(arise_chem_tests
    PRIVATE
        arise::chem
        arise::io
        arise_test_fixtures
        GTest::gtest_main)
arise_set_warnings(arise_chem_tests)

gtest_discover_tests(arise_chem_tests)
