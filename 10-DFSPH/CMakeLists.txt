cmake_minimum_required(VERSION 3.26)
project(DFSPH)
set(CMAKE_CXX_STANDARD 17)

include(../FindHoudini.cmake)

set(TARGET_NAME SIM_DFSPH)
add_library(
        ${TARGET_NAME}
        SHARED
        hdk/Entrance.cpp
        hdk/SIM_DFSPH_Particles.cpp
        hdk/SIM_DFSPH_Particles.h
        hdk/GAS_DFSPH_Solver.cpp
        hdk/GAS_DFSPH_Solver.h
        hdk/utils.h
)
target_link_libraries(
        ${TARGET_NAME}
        PUBLIC
        Houdini
)
houdini_configure_target(${TARGET_NAME})
