cmake_minimum_required(VERSION 3.18)

project(hypercurve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#############################################
##### 	Build libsndfile
#############################################
set(ENABLE_EXTERNAL_LIBS OFF CACHE BOOL "Build libsndfile external libs")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build Shared libs")
add_subdirectory(libsndfile)
target_compile_options(sndfile PRIVATE -fPIC)

#############################################
##### 	Build Hypercurve library
#############################################
add_library(hypercurve
   src/core.h
   src/curve_lib.h
   src/utilities.h
   src/cubic_spline.h

   # Include asciiplot library
   src/asciiplot/asciiplotter.h
   src/asciiplot/asciiplotter.cpp
   )
set_target_properties(hypercurve PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(hypercurve PRIVATE -fPIC)

#############################################
##### 	Build Test
#############################################
add_executable(hypercurve_test
   tests/test.cpp
   )

target_include_directories(hypercurve_test PUBLIC
   src/
   src/asciiplot
   )

target_link_libraries(hypercurve_test
   hypercurve
   sndfile
   )

#############################################
##### 	Build Lua Module
#############################################

if(NOT DEFINED LUA_INCLUDE_DIR OR NOT DEFINED LUA_LIBRARIES)
   find_package(Lua REQUIRED)
else()
   set(Lua_FOUND TRUE)
endif()

if(NOT ${Lua_FOUND})
   message(FATAL_ERROR "Lua not found")
endif()

add_library(lua_hypercurve SHARED
   lua_module/lua_hypercurve.cpp
   )
target_include_directories(lua_hypercurve PUBLIC
   ${LUA_INCLUDE_DIR}
   lua_module/lua-compat-5.3/c-api/
   )
target_link_libraries(lua_hypercurve hypercurve sndfile ${LUA_LIBRARIES})
set_target_properties(lua_hypercurve PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(lua_hypercurve PRIVATE -fPIC)

## Lua C Test

add_executable(lua_test lua_module/lua_test.cpp)
target_link_libraries(lua_test lua_hypercurve)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lua_module/test.lua
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
