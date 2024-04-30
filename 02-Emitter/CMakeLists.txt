cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(Emitter LANGUAGES CXX)
include(../FindHoudini.cmake)

add_library(
        SIM_Emitter
        SHARED
        hdk/Entrance.cpp
        hdk/GAS_Volume_Emitter.cpp
        hdk/GAS_Volume_Emitter.h
)
target_link_libraries(
        SIM_Emitter
        PUBLIC
        Houdini
)
target_include_directories(SIM_Emitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

houdini_get_default_install_dir(HOUDINI_INSTALL_DIR)
file(COPY DOP_Volume_Emitter.png DESTINATION ${HOUDINI_INSTALL_DIR}/config/icons)
houdini_configure_target(SIM_Emitter)
