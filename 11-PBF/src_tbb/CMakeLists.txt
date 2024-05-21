cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(PBF_TBB)

add_library(PBF_TBB STATIC PBF.cpp PBF.h)
target_link_libraries(PBF_TBB PUBLIC Houdini)


add_executable(PBF_TBB_test PBF.cpp PBF.h)
target_link_libraries(PBF_TBB_test PUBLIC Houdini)
target_compile_definitions(PBF_TBB_test PUBLIC TEST_PBF)
