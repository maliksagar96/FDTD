cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(FDTD_1D)

set(target_name fdtd_1d)

set(SOURCES
    01_fdtd_1d.cpp
    02_fdtd_graded_pml.cpp
    03_fdtd_graded_pml.cpp
    04_fdtd_CMPL.cpp
)

add_executable(${target_name} ${SOURCES})
target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})



