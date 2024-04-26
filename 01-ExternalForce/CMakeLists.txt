cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(ExternalForce LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_ExternalForce
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_Gravity_Force.cpp
        hdk/GAS_Gravity_Force.h
)
target_link_libraries(
        SIM_ExternalForce
        PUBLIC
        Houdini
)
target_include_directories(SIM_ExternalForce PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
houdini_configure_target(SIM_ExternalForce)
