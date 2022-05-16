cmake_minimum_required(VERSION 3.19)

execute_process(COMMAND quom ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_faust.h ${CMAKE_BINARY_DIR}/faust_lib/hypercurve_faust.h)
project(faust_library)

add_library(hypercurve_header INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve_faust.h ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.lib ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.dsp)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.dsp
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hypercurve.lib
   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

execute_process(COMMAND faust ${CMAKE_CURRENT_BINARY_DIR}/hypercurve.dsp -o ${CMAKE_CURRENT_BINARY_DIR}/hh.cpp)
execute_process(COMMAND faust2alsa ${CMAKE_CURRENT_BINARY_DIR}/hypercurve.dsp)
