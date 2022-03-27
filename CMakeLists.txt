cmake_minimum_required(VERSION 3.5)

project(diocles LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(diocles main.cpp)

target_link_libraries(diocles sndfile)
