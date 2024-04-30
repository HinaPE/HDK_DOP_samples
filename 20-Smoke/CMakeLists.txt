cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(Smoke LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_Smoke
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_Smoke_Solver.cpp
        hdk/GAS_Smoke_Solver.h
)
target_link_libraries(
        SIM_Smoke
        PUBLIC
        Houdini
)
target_include_directories(SIM_Smoke PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
file(COPY DOP_Smoke_Solver.png DESTINATION ${HOUDINI_INSTALL_DIR}/config/icons)
houdini_configure_target(SIM_Smoke)
