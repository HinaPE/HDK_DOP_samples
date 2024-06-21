cmake_minimum_required(VERSION 3.26)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(StableFluids)
include(../FindHoudini.cmake)

add_library(
        SIM_StableFluids
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_StableFluids.cpp
        hdk/GAS_StableFluids.h
        hdk/poisson.cpp
        hdk/poisson.h
        hdk/diffusion.cpp
        hdk/diffusion.h
)
target_link_libraries(
        SIM_StableFluids
        PUBLIC
        Houdini
)
target_include_directories(SIM_StableFluids PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
configure_file(../icon.png ${HOUDINI_INSTALL_DIR}/config/icons/DOP_StableFluids.png COPYONLY)
houdini_configure_target(SIM_StableFluids)
