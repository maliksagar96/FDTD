cmake_minimum_required(VERSION 3.10)
project(Mesh2D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(VTK REQUIRED)

include_directories(/usr/include/opencascade include)

set(SOURCES
    src/main.cpp
    src/Mesh2D.cpp
    src/vtkWriter.cpp
)

add_executable(read_step ${SOURCES})

set(OCC_LIBS
    TKernel
    TKMath
    TKBRep
    TKSTEP
    TKSTEPAttr
    TKSTEP209
    TKSTEPBase
    TKIGES
    TKXSBase
    TKGeomBase
    TKGeomAlgo
    TKTopAlgo
)

target_link_libraries(read_step PRIVATE
    ${VTK_LIBRARIES}
    ${OCC_LIBS}
)
