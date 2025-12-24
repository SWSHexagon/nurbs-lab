cmake_minimum_required(VERSION 3.16)
project(nurbs_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32 AND DEFINED WINVER_VALUE)
    add_compile_definitions(_WIN32_WINNT=${WINVER_VALUE})
endif()

# vcpkg makes these packages available automatically
find_package(Boost REQUIRED COMPONENTS thread)

option(BUILD_TESTS "Build unit tests" ON)

add_library(nurbs
    src/surface_curvature.hpp
    src/surface_differential.hpp
    src/surface_frame.hpp
    src/closest_point_result.hpp
    src/bspline_basis.cpp src/bspline_basis.hpp
    src/bspline_curve.cpp src/bspline_curve.hpp
    src/bspline_surface.cpp src/bspline_surface.hpp
    src/nurbs_surface.cpp src/nurbs_surface.hpp
    src/math_utils.cpp src/math_utils.hpp
    src/bspline_surface_builder.hpp
    src/bspline_surface_sphere.cpp
    src/bspline_surface_cylinder.cpp
    src/bspline_surface_plane.cpp
    src/bspline_surface_bowl.cpp
    src/bspline_curve_builder.hpp
    src/bspline_curve_sample.cpp
    src/bspline_curve_stress.cpp
    src/parametric_curve.hpp
    src/parametric_surface.hpp
    src/parametric_surface.cpp
    src/parametric_curve.cpp
)

target_include_directories(nurbs PUBLIC src)

target_link_libraries(nurbs PRIVATE Boost::thread)

add_executable(nurbs_demo src/main.cpp)
target_link_libraries(nurbs_demo PRIVATE nurbs)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
