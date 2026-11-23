# Apps: the thin layers around the engine. apps/wasm exists only for Emscripten builds.
if(EMSCRIPTEN)
    add_subdirectory(wasm)
endif()
