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


if(DEFINED BUILD_ALL)
   set(BUILD_CSOUND_OPCODE TRUE)
   set(BUILD_LUA_MODULE TRUE)
endif()


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

if(${BUILD_LUA_MODULE})
   message(STATUS "Building Lua module")
   add_subdirectory(lua_module)
   target_link_libraries(lua_hypercurve PUBLIC hypercurve_static sndfile ${LUA_LIBRARIES})
endif()

#############################################
##### 	Build Csound Opcode
#############################################

if(${BUILD_CSOUND_OPCODE})
   message(STATUS "Building Csound opcode")
   add_subdirectory(csound_opcode)
   target_link_libraries(csound_hypercurve PUBLIC hypercurve_static)
endif()

