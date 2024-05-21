project(raytracer)
cmake_minimum_required(VERSION 3.17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED 17)

find_package(OpenMP)

add_executable(raytracer raytracer.cpp)
add_library(mesh STATIC mesh.cpp)
target_link_libraries(raytracer PRIVATE OpenMP::OpenMP_CXX PRIVATE mesh)

target_include_directories(raytracer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/stb)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")





