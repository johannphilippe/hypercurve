cmake_minimum_required(VERSION 2.8.12)
project(Csound_Hypercurve)

# Release or Debug
set(CMAKE_BUILD_TYPE "Release")

# switch RPATH use off
set(CMAKE_MACOSX_RPATH 0)
# force make to print the command lines
set(CMAKE_VERBOSE_MAKEFILE on)
# path to Csound cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
"${CMAKE_SOURCE_DIR}/../cmake/Modules/")
# set C++ compilation flags
set(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS}  "-std=c++17")

# Csound opcode build
find_package(Csound REQUIRED)

message(STATUS "CSOUND INCLUDE DIRS : ${CSOUND_INCLUDE_DIRS}")
message(STATUS "CSOUND LIBRARIES : ${CSOUND_LIBRARIES}")

function(make_plugin pluginlib srcs)
  if(APPLE)
        add_library(${pluginlib} SHARED ${srcs})
    else()
        add_library(${pluginlib} MODULE ${srcs})
  endif()
  target_include_directories(${pluginlib} PRIVATE ${CSOUND_INCLUDE_DIRS})
  target_compile_options(${pluginlib} PRIVATE -O3)
endfunction(make_plugin)

message(STATUS "Building Opcodes")
make_plugin(csound_hypercurve csound_hypercurve.cpp)

if(BUILD_TESTS)
      add_executable(csound_test test_csound.cpp)
      target_link_libraries(csound_test PUBLIC ${CSOUND_LIBRARIES})
      target_include_directories(csound_test PUBLIC ${CSOUND_INCLUDE_DIRS})
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test.csd
         DESTINATION ${CMAKE_BINARY_DIR}/bin)
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypertest.csd
         DESTINATION ${CMAKE_BINARY_DIR}/bin)
endif()
# Installation
#[[
if(NOT DEFINED CSOUND_INSTALL_PATH)
    set(CSOUND_INSTALL_PATH ${CSOUND_PLUGINS_DIR})
endif()
message(STATUS "Csound Install path : ${CSOUND_INSTALL_PATH}")
install(TARGETS csound_hypercurve
    RUNTIME DESTINATION "${CSOUND_INSTALL_PATH}"
    LIBRARY DESTINATION "${CSOUND_INSTALL_PATH}"
    )
]]
