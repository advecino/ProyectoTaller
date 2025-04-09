cmake_minimum_required(VERSION 3.28)
project(ProyectoTaller)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(ProyectoTaller
        src/Matrix.cpp
        include/Matrix.h
        src/Mjday.cpp
        include/Mjday.h
        src/Frac.cpp
        include/Frac.h
        src/gmst.cpp
        include/gmst.h
        EKF_Test.cpp
        src/R_x.cpp
        include/R_x.h
        include/Sat_const.h
        EKF_GEOS3.cpp
)
