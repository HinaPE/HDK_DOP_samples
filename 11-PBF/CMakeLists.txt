cmake_minimum_required(VERSION 3.26)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O3 /GL /fp:fast /openmp /arch:AVX2 /Qvec") # MSVC
if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
endif ()
project(PBF LANGUAGES CXX CUDA)
include(../FindHoudini.cmake)

add_subdirectory(src_tbb)

add_library(
        SIM_PBF
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_PBF_Solver.cpp
        hdk/GAS_PBF_Solver.h
)
target_link_libraries(
        SIM_PBF
        PUBLIC
        Houdini
        PBF_TBB
)
target_include_directories(SIM_PBF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
configure_file(../icon.png ${HOUDINI_INSTALL_DIR}/config/icons/DOP_PBF_Solver.png COPYONLY)
houdini_configure_target(SIM_PBF)
