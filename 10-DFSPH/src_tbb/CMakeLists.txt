cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(DFSPH_TBB)

add_library(DFSPH_TBB STATIC DFSPH.cpp DFSPH.h)
target_link_libraries(DFSPH_TBB PUBLIC Houdini)


add_executable(DFSPH_TBB_test DFSPH.cpp DFSPH.h)
target_link_libraries(DFSPH_TBB_test PUBLIC Houdini)
target_compile_definitions(DFSPH_TBB_test PUBLIC TEST_DFSPH)
