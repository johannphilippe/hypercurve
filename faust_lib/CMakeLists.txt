cmake_minimum_required(VERSION 3.16)

execute_process(COMMAND quom
   ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_faust.h
   ${CMAKE_CURRENT_SOURCE_DIR}/faust_hypercurve_lib/hypercurve_faust.h)

project(faust_library)

#add_library(hypercurve_header INTERFACE
#   ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_faust.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.lib
#   ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_test.dsp)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_test.dsp
   DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/faust_hypercurve_lib)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.lib
   DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/faust_hypercurve_lib)

# Copy files to binary dir for tests
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_test.dsp
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.lib
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/faust_hypercurve_lib/hypercurve_faust.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Commands to test with faust
#execute_process(COMMAND faust ${CMAKE_CURRENT_BINARY_DIR}/hypercurve.dsp -o ${CMAKE_CURRENT_BINARY_DIR}/hypercurve_faust.cpp)
#execute_process(COMMAND faust2alsa ${CMAKE_CURRENT_BINARY_DIR}/hypercurve.dsp)
