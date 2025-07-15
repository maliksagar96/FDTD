cmake_minimum_required(VERSION 3.20)

project(FDTD_1D)
set(CMAKE_CXX_STANDARD 17)

set(target_name FDTD_1D)

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED jsoncpp)

add_executable(${target_name}
    src/main.cpp
    src/FDTD_1D.cpp
)

include_directories(include ${JSONCPP_INCLUDE_DIRS})
target_link_libraries(${target_name} PRIVATE ${JSONCPP_LIBRARIES})
