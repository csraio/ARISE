add_library(arise_pipeline
    json.cpp
    pipeline.cpp)
add_library(arise::pipeline ALIAS arise_pipeline)

target_link_libraries(arise_pipeline PUBLIC arise::filter arise::io)
arise_set_warnings(arise_pipeline)
