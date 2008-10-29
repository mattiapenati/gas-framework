# Common
set(SOURCES Common/Exception.cpp Common/Iterable.cpp Common/Program.cpp)

# LinearAlgebra
set(SOURCES ${SOURCES} LinearAlgebra/Vector.cpp LinearAlgebra/Matrix.cpp LinearAlgebra/Solver.cpp LinearAlgebra/SparseMatrix.cpp)


add_library(Gas SHARED ${SOURCES})
