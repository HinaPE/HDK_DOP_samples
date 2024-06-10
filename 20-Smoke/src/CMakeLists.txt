cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(SmokeSrc)

add_library(SmokeSrc STATIC stable_fluids.cpp stable_fluids.h)
target_link_libraries(SmokeSrc PUBLIC Houdini)

