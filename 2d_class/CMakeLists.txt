cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(fdtd_2d)

set(target_name fdtd_2d)


find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

set(sources
    src/fdtd_2d.cpp
    src/main.cpp
    src/publish.cpp
)

include_directories(${JSONCPP_INCLUDE_DIRS} include)
link_directories(${JSONCPP_LIBRARY_DIRS})
add_executable(${target_name} ${sources})

target_link_libraries(${target_name} ${JSONCPP_LIBRARIES})
