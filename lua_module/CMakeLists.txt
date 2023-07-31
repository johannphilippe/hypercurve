cmake_minimum_required(VERSION 3.18)
project(Lua_Hypercurve LANGUAGES CXX)

if(NOT DEFINED LUA_INCLUDE_DIR OR NOT DEFINED LUA_LIBRARIES)
   find_package(Lua REQUIRED)
else()
   set(Lua_FOUND TRUE)
endif()

message(STATUS ${LUA_VERSION_STRING})
message(STATUS "${LUA_INCLUDE_DIR}")
message(STATUS "${LUA_LIBRARIES}")

if(NOT ${Lua_FOUND})
   message(FATAL_ERROR "Lua not found")
endif()

message(STATUS "Shared exports generated")

add_library(lua_hypercurve SHARED
   lua_hypercurve.cpp
   lua-compat-5.3/c-api/compat-5.3.h
   )

include (GenerateExportHeader)
GENERATE_EXPORT_HEADER(lua_hypercurve
    BASE_NAME lua_hypercurve
    EXPORT_MACRO_NAME SHARED_EXPORTS
    EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/SHARED_EXPORTS.h
    STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)

target_include_directories(lua_hypercurve PUBLIC
   ${LUA_INCLUDE_DIR}
   lua-compat-5.3/c-api/
   ${CMAKE_BINARY_DIR}/exports
   )

target_link_libraries(lua_hypercurve PUBLIC ${LUA_LIBRARIES} hypercurve_static)
target_compile_options(lua_hypercurve PUBLIC -fPIC)

if(WIN32)
   set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
   if(${CMAKE_BUILD_TYPE} MATCHES "Release")
       message(STATUS "Build Release")
       set(WIN_COMPILE_FLAGS "/c /nologo /MD /TP /Dx64")
   else()
       message(STATUS "Build debug")
       set(WIN_COMPILE_FLAGS "/c /nologo /MDd /TP /Dx64")
   endif()
   set(WIN_LINK_FLAGS "/nologo /dll")
   set_target_properties(lua_hypercurve
      PROPERTIES
      LINKER_LANGUAGE CXX
      COMPILE_FLAGS ${WIN_COMPILE_FLAGS}
      LINK_FLAGS ${WIN_LINK_FLAGS}
      PREFIX ""
      OUTPUT_NAME "liblua_hypercurve"
      )
else()

   set_target_properties(lua_hypercurve
      PROPERTIES
      LINKER_LANGUAGE CXX
      PREFIX "lib"
      OUTPUT_NAME "lua_hypercurve"
      )
endif()

## Lua CPP Test
if(BUILD_TESTS)(
      add_executable(lua_test lua_test.cpp)
      target_link_libraries(lua_test ${LUA_LIBRARIES})
      target_include_directories(lua_test PUBLIC ${LUA_INCLUDE_DIR})
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test.lua
         DESTINATION ${CMAKE_BINARY_DIR}/bin)
)
if(WIN32 AND NOT BUILD_LUAJIT)
    string(REPLACE ".lib" ".dll" LUA_DLL ${LUA_LIBRARIES})
    message(STATUS ${LUA_DLL})
    file(COPY ${LUA_DLL} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endif()
