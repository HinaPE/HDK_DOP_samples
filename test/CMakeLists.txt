cmake_minimum_required(VERSION 3.26)
project(TEST)
set(CMAKE_CXX_STANDARD 20)

set(TARGET_NAME MY_TEST)
add_executable(${TARGET_NAME} main.cpp)
target_include_directories(${TARGET_NAME} PRIVATE ../3rd_party/base/)
