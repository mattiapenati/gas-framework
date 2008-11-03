# Common
set(SOURCES
	Common/Exception.cpp
	Common/Iterable.cpp
	Common/Program.cpp
)

# LinearAlgebra
set(SOURCES ${SOURCES}
	LinearAlgebra/Matrix.cpp
	LinearAlgebra/Solver.cpp
	LinearAlgebra/SparseMatrix.cpp
)

# Library
add_library(Gas ${SOURCES})
