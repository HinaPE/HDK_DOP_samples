cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(VolumeRenderer LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_VolumeRenderer
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_Volume_Renderer.cpp
        hdk/GAS_Volume_Renderer.h
)
target_link_libraries(
        SIM_VolumeRenderer
        PUBLIC
        Houdini
)
target_include_directories(SIM_VolumeRenderer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
configure_file(../icon.png ${HOUDINI_INSTALL_DIR}/config/icons/DOP_Volume_Renderer.png COPYONLY)
houdini_configure_target(SIM_VolumeRenderer)
