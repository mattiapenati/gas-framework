include_directories(
	${GAS_SOURCE_DIR}/Gas/SparseLib++/
	${GAS_SOURCE_DIR}/Gas/SparseLib++/mv/
)
add_definitions(-DCOMPLEX=std::complex<double>)


add_library(sparse STATIC
	compcol_double.cc
	comprow_double.cc
	coord_double.cc
	iotext_double.cc
	diagpre_double.cc
	icpre_double.cc
	ilupre_double.cc
	qsort_double.cc
	qsort_int.cc
	iohb_double.cc
	iohb.c
)

add_subdirectory(spblas)
add_subdirectory(mv)
