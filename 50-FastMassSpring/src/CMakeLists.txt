cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
project(FastMassSpringImpl)

add_library(FastMassSpringImpl STATIC FastMassSpring.cpp FastMassSpring.h)
target_link_libraries(FastMassSpringImpl PUBLIC Houdini)


add_executable(FastMassSpringImpl_test FastMassSpring.cpp FastMassSpring.h)
target_link_libraries(FastMassSpringImpl_test PUBLIC Houdini)
target_compile_definitions(FastMassSpringImpl_test PUBLIC TEST_FastMassSpring)
