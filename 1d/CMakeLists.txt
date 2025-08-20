cmake_minimum_required(VERSION 3.10)
project(FDTD_1D CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sources
set(SOURCES
    src/main.cpp
    src/fdtd_1d.cpp
)

# Include directories
include_directories(include)

# Find JsonCpp (system package)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

if (JSONCPP_FOUND)
    include_directories(${JSONCPP_INCLUDE_DIRS})
    link_directories(${JSONCPP_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "JsonCpp not found! Please install libjsoncpp-dev")
endif()

# Executable
add_executable(fdtd ${SOURCES})

# Link libraries
target_link_libraries(fdtd ${JSONCPP_LIBRARIES})
