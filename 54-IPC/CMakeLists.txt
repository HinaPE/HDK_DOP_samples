cmake_minimum_required(VERSION 3.26)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /fp:fast /openmp /arch:AVX2") # MSVC
project(IPC LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_IPC
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_IPC_Solver.cpp
        hdk/GAS_IPC_Solver.h
)
target_link_libraries(
        SIM_IPC
        PUBLIC
        Houdini
)

target_include_directories(SIM_IPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
configure_file(../icon.png ${HOUDINI_INSTALL_DIR}/config/icons/DOP_IPC_Solver.png COPYONLY)
houdini_configure_target(SIM_IPC)
