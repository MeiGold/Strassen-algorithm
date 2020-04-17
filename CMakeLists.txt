cmake_minimum_required(VERSION 3.16)
project(StrassenAlgorithm)

set(CMAKE_CXX_STANDARD 14)

add_executable(StrassenAlgorithm main.cpp)
target_link_libraries(StrassenAlgorithm pthread)