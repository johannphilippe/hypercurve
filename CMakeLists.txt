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
   )
set_target_properties(hypercurve PROPERTIES LINKER_LANGUAGE CXX)

#############################################
##### 	Build Test
#############################################
add_executable(hypercurve_test
   tests/test.cpp
   )

target_include_directories(hypercurve_test PUBLIC
   src/core.h
   src/curve_lib.h
   )

target_link_libraries(hypercurve_test
   sndfile
   )

#############################################
##### 	Build Lua Module
#############################################
find_package(Lua REQUIRED)

if(NOT ${Lua_FOUND})
   message(FATAL_ERROR "Lua not found")
endif()


add_library(lua_hypercurve SHARED lua_module/lua_hypercurve.cpp)
target_include_directories(lua_hypercurve PUBLIC ${LUA_INCLUDE_DIR})
target_link_libraries(lua_hypercurve hypercurve sndfile ${LUA_LIBRARIES})
set_target_properties(lua_hypercurve PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(lua_hypercurve PRIVATE -fPIC)
