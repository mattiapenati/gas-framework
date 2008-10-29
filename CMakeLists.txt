cmake_minimum_required(VERSION 2.6)
project(GAS)
include_directories(${GAS_SOURCE_DIR})
link_directories (${GAS_BINARY_DIR}/Gas)

# add_subdirectory(Documentation)
add_subdirectory(Gas)
# add_subdirectory(Main)
# add_subdirectory(Test)
