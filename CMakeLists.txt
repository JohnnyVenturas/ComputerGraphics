cmake_minimum_required(VERSION 3.17)
project(simulator)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED 17)

find_package(OpenMP)

add_executable(simulator main.cpp)
add_library(mesh STATIC mesh.cpp)
add_library(Ray STATIC Ray.cpp)
add_library(Scene STATIC Scene.cpp)
add_library(Vector STATIC Vector.cpp)
add_library(raytracer STATIC raytracer.cpp)
add_library(fluids STATIC fluids.cpp)
add_library(lbfgs STATIC lbfgs.c)
target_include_directories(lbfgs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
 #add_library(Camera INTERFACE camera.h)
 #add_library(Intersect INTERFACE Intersect.h)

target_link_libraries(simulator PRIVATE OpenMP::OpenMP_CXX PRIVATE mesh Ray Scene Vector raytracer fluids)
#target_link_libraries(polygon PRIVATE Vector)
target_link_libraries(fluids PRIVATE Vector lbfgs OpenMP::OpenMP_CXX)

target_include_directories(fluids PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/stb ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(simulator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/stb ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(raytracer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/stb ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wredundant-move -Wpessimizing-move")





