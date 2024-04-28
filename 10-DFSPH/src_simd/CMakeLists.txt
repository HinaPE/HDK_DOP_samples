cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(DFSPH_SIMD)

add_library(DFSPH_SIMD STATIC DFSPH.cpp DFSPH.h)
target_link_libraries(DFSPH_SIMD PUBLIC Houdini)


add_executable(DFSPH_SIMD_test DFSPH.cpp DFSPH.h)
target_link_libraries(DFSPH_SIMD_test PUBLIC Houdini)
target_compile_definitions(DFSPH_SIMD_test PUBLIC TEST_DFSPH)
