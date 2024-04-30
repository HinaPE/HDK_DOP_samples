cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(CHOP LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_CHOP
        SHARED
        hdk/Entrance.cpp
        hdk/COP2_Stochastic_Tomography.cpp
        hdk/COP2_Stochastic_Tomography.cpp
)
target_link_libraries(
        SIM_CHOP
        PUBLIC
        Houdini
)
target_include_directories(SIM_CHOP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
houdini_configure_target(SIM_CHOP)
