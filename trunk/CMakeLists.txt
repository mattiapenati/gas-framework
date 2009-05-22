cmake_minimum_required(VERSION 2.6)
project(GAS CXX)

include_directories(${GAS_SOURCE_DIR})
link_directories(${GAS_BINARY_DIR}/Gas)

# This only for SparseLib++
include_directories(
        ${GAS_SOURCE_DIR}/Gas/SparseLib++/
        ${GAS_SOURCE_DIR}/Gas/SparseLib++/mv/
		/opt/local/include
)
add_definitions(-DCOMPLEX=std::complex<double>)
set(SPARSELIB mv sparse spblas)

# Subdirectory
add_subdirectory(Documentation)
add_subdirectory(Gas)
add_subdirectory(Gas/SparseLib++)
add_subdirectory(Main)
add_subdirectory(Test)
