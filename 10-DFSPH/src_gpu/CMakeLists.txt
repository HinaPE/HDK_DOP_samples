cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)
if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
endif ()
project(DFSPH_GPU LANGUAGES CXX CUDA)


file(GLOB_RECURSE NEIGHBOR_SEARCH_SOURCE neighbor/*.h neighbor/*.cpp neighbor/*.cu)
add_library(NeighborSearch_GPU STATIC ${NEIGHBOR_SEARCH_SOURCE})
target_include_directories(NeighborSearch_GPU PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include(GenerateExportHeader)
generate_export_header(NeighborSearch_GPU BASE_NAME cuNSearch EXPORT_MACRO_NAME cuNSearch_EXPORT EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/cuNSearch/cuNSearch_export.h)
target_include_directories(NeighborSearch_GPU PUBLIC neighbor/include neighbor/Utils)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/cuNSearch/cuNSearch_export.h DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/neighbor/include)


add_library(DFSPH_GPU STATIC DFSPH.cu DFSPH.cuh)
target_link_libraries(DFSPH_GPU PUBLIC NeighborSearch_GPU)
target_compile_options(DFSPH_GPU PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
target_compile_options(DFSPH_GPU PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)


add_executable(DFSPH_GPU_test DFSPH.cu DFSPH.cuh)
target_link_libraries(DFSPH_GPU_test PUBLIC NeighborSearch_GPU)
target_compile_options(DFSPH_GPU_test PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
target_compile_options(DFSPH_GPU_test PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
target_compile_definitions(DFSPH_GPU_test PUBLIC TEST_DFSPH)
