cmake_minimum_required(VERSION 3.19)

project(hypercurve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

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
set(HYPERCURVE_SOURCES
   src/core.h
   src/curve_lib.h
   src/utilities.h
   src/cubic_spline.h
   src/modulator_lib.h
   src/hypercurve.h

   # Include asciiplot library
   src/asciiplot/asciiplotter.h
   src/asciiplot/asciiplotter.cpp
   )

add_library(hypercurve SHARED
   ${HYPERCURVE_SOURCES}
   )
set_target_properties(hypercurve PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(hypercurve PRIVATE -fPIC)

add_library(hypercurve_static STATIC
   ${HYPERCURVE_SOURCES}
   )
set_target_properties(hypercurve_static PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(hypercurve_static PRIVATE -fPIC)


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

message(STATUS ${LUA_VERSION_STRING})
message(STATUS ${LUA_INCLUDE_DIR})
message(STATUS ${LUA_LIBRARIES})

if(NOT ${Lua_FOUND})
   message(FATAL_ERROR "Lua not found")
endif()

add_library(lua_hypercurve SHARED
   lua_module/lua_hypercurve.cpp
   lua_module/lua-compat-5.3/c-api/compat-5.3.h
   )
target_include_directories(lua_hypercurve PUBLIC
   ${LUA_INCLUDE_DIR}
   lua_module/lua-compat-5.3/c-api/
   )

target_link_libraries(lua_hypercurve PUBLIC hypercurve sndfile ${LUA_LIBRARIES})
set_target_properties(lua_hypercurve PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(lua_hypercurve PUBLIC -fPIC)

## Lua C Test

add_executable(lua_test lua_module/lua_test.cpp)
target_link_libraries(lua_test lua_hypercurve)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lua_module/test.lua
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


#############################################
##### 	Build Csound Opcode
#############################################

if(NOT DEFINED CSOUND_INCLUDE_DIRS OR NOT DEFINED CSOUND_LIBRARIES)
   find_package(Csound REQUIRED)
else()
   set(CSOUND_FOUND TRUE)
endif()

message(STATUS ${CSOUND_INCLUDE_DIRS})
message(STATUS ${CSOUND_LIBRARIES})

if(CSOUND_FOUND)
   add_library(csound_hypercurve SHARED csound_opcode/csound_hypercurve.cpp)
   target_link_libraries(csound_hypercurve PUBLIC hypercurve_static)
   target_include_directories(csound_hypercurve PRIVATE ${CSOUND_INCLUDE_DIRS})
   target_compile_options(csound_hypercurve PUBLIC -fPIC)
endif()

## Csound Cpp Test

add_executable(csound_test csound_opcode/test_csound.cpp)
target_link_libraries(csound_test PUBLIC ${CSOUND_LIBRARIES})
target_include_directories(csound_test PUBLIC ${CSOUND_INCLUDE_DIRS})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/csound_opcode/test.csd
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/csound_opcode/other_test.csd
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
