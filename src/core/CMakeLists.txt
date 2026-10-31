add_library(arise_core
    error.cpp
    formula.cpp)
add_library(arise::core ALIAS arise_core)

target_include_directories(arise_core
    PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_compile_features(arise_core PUBLIC cxx_std_20)
arise_set_warnings(arise_core)
