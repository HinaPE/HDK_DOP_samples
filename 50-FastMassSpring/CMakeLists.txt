cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(FastMassSpring LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_FastMassSpring
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_FMS_Solver.cpp
        hdk/GAS_FMS_Solver.h
)
target_link_libraries(
        SIM_FastMassSpring
        PUBLIC
        Houdini
)
target_include_directories(SIM_FastMassSpring PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
file(COPY DOP_FMS_Solver.png DESTINATION ${HOUDINI_INSTALL_DIR}/config/icons)
houdini_configure_target(SIM_FastMassSpring)
