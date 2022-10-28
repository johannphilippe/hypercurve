cmake_minimum_required(VERSION 3.19)

project(hypercurve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "Release")
    #add_compile_definitions($<$<CONFIG:Release>:_ITERATOR_DEBUG_LEVEL=0>)
else()
    message(STATUS "Debug")
    #add_compile_definitions($<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=1>)

endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

#############################################
##### 	Options
#############################################

option(BUILD_CSOUND_OPCODE "Build Csound opcodes" ON)
option(BUILD_LUA_MODULE "Build Lua module" ON)
option(BUILD_LUAJIT "Build Luajit" OFF)
option(BUILD_FAUST_LIB	"Build Faust library" ON)
option(BUILD_ALL "Build all libraries" ON)
option(SSE, "SSE support for FPNG" ON)
option(MACOS_UNIVERSAL "Support for MacOS Universal build - both x86_64 and ARM64" ON)

if(MACOS_UNIVERSAL) 
	set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")
endif()

#############################################
##### 	Build libsndfile
#############################################
set(ENABLE_EXTERNAL_LIBS OFF CACHE BOOL "Build libsndfile external libs")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build Shared libs")
set(BUILD_EXAMPLES OFF CACHE BOOL "Not Building Examples for libsndfile")
set(BUILD_PROGRAMS OFF CACHE BOOL "Not Building Programs for libsndfile")
set(BUILD_TESTING OFF CACHE BOOL "Not Build testing for libsndfile")
#set(ENABLE_STATIC_RUNTIME ON CACHE BOOL "Enable static runtime")
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
   src/extra.h
   src/hypercurve.h

   # fpng library
   src/fpng/src/fpng.h
   src/fpng/src/fpng.cpp

   # asciiplot library
   src/asciiplot/asciiplotter.h
   src/asciiplot/asciiplotter.cpp
   )


add_library(hypercurve SHARED
   ${HYPERCURVE_SOURCES}
   )

set_target_properties(hypercurve PROPERTIES LINKER_LANGUAGE CXX)
# To disable SSE use the following and remove -msse4.1 -mpclmul
set(HYPERCURVE_OPTIONS -fPIC)
if(SSE)
   list(APPEND HYPERCURVE_OPTIONS -msse4.1 -mpclmul)
endif()

target_compile_options(hypercurve
   PRIVATE
   ${HYPERCURVE_OPTIONS}
   )

# Static version
add_library(hypercurve_static STATIC
   ${HYPERCURVE_SOURCES}
   )
set_target_properties(hypercurve_static PROPERTIES LINKER_LANGUAGE CXX)

# To disable SSE use the following and remove -msse4.1 -mpclmul
target_compile_definitions(hypercurve_static PUBLIC FPNG_NO_SSE)
target_compile_options(hypercurve_static
   PRIVATE
   ${HYPERCURVE_OPTIONS}
   )

if(NOT SSE)
   target_compile_definitions(hypercurve PUBLIC FPNG_NO_SSE)
   target_compile_definitions(hypercurve_static PUBLIC FPNG_NO_SSE)
endif()


#############################################
##### 	Build Test
#############################################

add_executable(hypercurve_test
   tests/test.cpp
   )

target_include_directories(hypercurve_test PUBLIC
   src/
   src/asciiplot
   src/fpng/src
   )

target_link_libraries(hypercurve_test
   hypercurve_static
   sndfile
   )

#############################################
##### 	Build Lua Module
#############################################

if(BUILD_LUA_MODULE OR BUILD_ALL)

   if(BUILD_LUAJIT)
      add_subdirectory(LuaJIT)
      target_compile_options(libluajit PRIVATE -fPIC)
      set(LUA_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/LuaJIT/src)
      set(LUA_LIBRARIES libluajit)
   endif()

   message(STATUS "Building Lua module ON")
   add_subdirectory(lua_module)
   #target_link_libraries(lua_hypercurve PUBLIC hypercurve_static sndfile)
endif()

#############################################
##### 	Build Csound Opcode
#############################################

if(BUILD_CSOUND_OPCODE OR BUILD_ALL)
   message(STATUS "Building Csound opcode ON")
   add_subdirectory(csound_opcode)
   target_link_libraries(csound_hypercurve PUBLIC hypercurve_static)
endif()


#############################################
##### 	Build Faust lib
#############################################

if(BUILD_FAUST_LIB OR BUILD_ALL)
   message(STATUS "Building Faust library ON")
   add_subdirectory(faust_lib)
endif()
