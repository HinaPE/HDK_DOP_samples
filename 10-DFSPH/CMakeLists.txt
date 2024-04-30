cmake_minimum_required(VERSION 3.26)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_CXX_STANDARD 17)
if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
endif ()
project(DFSPH LANGUAGES CXX CUDA)
include(../FindHoudini.cmake)

add_subdirectory(src_tbb)
add_subdirectory(src_simd)
add_subdirectory(src_gpu)

add_library(
        SIM_DFSPH
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_DFSPH_Solver.cpp
        hdk/GAS_DFSPH_Solver.h
)
target_link_libraries(
        SIM_DFSPH
        PUBLIC
        Houdini
        DFSPH_TBB
        DFSPH_SIMD
        DFSPH_GPU
)
target_include_directories(SIM_DFSPH PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
configure_file(../icon.png ${HOUDINI_INSTALL_DIR}/config/icons/DOP_DFSPH_Solver.png COPYONLY)
houdini_configure_target(SIM_DFSPH)
