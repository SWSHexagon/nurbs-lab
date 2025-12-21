cmake_minimum_required(VERSION 3.16)
project(nurbs_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build unit tests" ON)

add_library(nurbs
    src/bspline_basis.cpp
    src/bspline_curve.cpp
    src/bspline_surface.cpp
    src/nurbs_surface.cpp
    src/math_utils.cpp
    src/test_closest_point.cpp
)

target_include_directories(nurbs PUBLIC src)

add_executable(nurbs_demo src/main.cpp)
target_link_libraries(nurbs_demo PRIVATE nurbs)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
