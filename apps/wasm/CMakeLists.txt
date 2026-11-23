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

# The web client as one folder to serve: the page (public/, without its tests) and the
# engine next to it, at build/<preset>/site. Serve that folder over HTTP (ES modules do not
# load from file://), e.g. python3 -m http.server --directory build/web/site 8000.
set(ARISE_SITE_DIR "${PROJECT_BINARY_DIR}/site")
add_custom_target(arise_site ALL
    COMMAND "${CMAKE_COMMAND}" -E rm -rf "${ARISE_SITE_DIR}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${ARISE_SITE_DIR}/engine"
    COMMAND "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/public/index.html" "${ARISE_SITE_DIR}/"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory "${PROJECT_SOURCE_DIR}/public/css" "${ARISE_SITE_DIR}/css"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory "${PROJECT_SOURCE_DIR}/public/js" "${ARISE_SITE_DIR}/js"
    COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:arise_engine>"
            "$<TARGET_FILE_DIR:arise_engine>/arise_engine.wasm" "${ARISE_SITE_DIR}/engine/"
    DEPENDS arise_engine
    COMMENT "Assembling the web client in ${ARISE_SITE_DIR}"
    VERBATIM)
