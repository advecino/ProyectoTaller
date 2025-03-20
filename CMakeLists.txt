cmake_minimum_required(VERSION 3.28)
project(ProyectoTaller)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(ProyectoTaller
        src/Matrix.cpp
        include/Matrix.h)
