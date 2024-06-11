cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(SmokeSrc)

add_library(SmokeSrc STATIC pressure.cpp pressure.h diffuse.cpp diffuse.h advect.cpp advect.h source.cpp source.h)
target_link_libraries(SmokeSrc PUBLIC Houdini)
