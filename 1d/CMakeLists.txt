cmake_minimum_required(VERSION 3.10)

project(fdtd_1d_simulation)

add_executable(fdtd_1d_simulation src/main.cpp src/FDTD_1D.cpp)

target_include_directories(fdtd_1d_simulation PRIVATE include)