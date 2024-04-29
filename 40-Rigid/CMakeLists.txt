cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(Rigid LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_Rigid
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_Rigid_Apply_Impulse.cpp
        hdk/GAS_Rigid_Apply_Impulse.h
)
target_link_libraries(
        SIM_Rigid
        PUBLIC
        Houdini
)
target_include_directories(SIM_Rigid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
houdini_configure_target(SIM_Rigid)
