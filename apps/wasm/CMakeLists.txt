# The engine for the web client: an ES module (arise_engine.js + arise_engine.wasm) whose
# functions take plain JavaScript values and return JSON text (see bindings.cpp).
add_executable(arise_engine bindings.cpp)
target_link_libraries(arise_engine PRIVATE arise::pipeline)
arise_set_warnings(arise_engine)
target_link_options(arise_engine PRIVATE
    -lembind
    -sMODULARIZE=1
    -sEXPORT_ES6=1
    -sEXPORT_NAME=createAriseEngine
    -sENVIRONMENT=web,node
    -sALLOW_MEMORY_GROWTH=1)
set_target_properties(arise_engine PROPERTIES SUFFIX ".js")
