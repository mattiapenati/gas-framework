FIND_PATH(
  EIGEN_INCLUDE_DIRS Eigen/Array
  /usr/include/eigen3
  /usr/local/include/eigen3
  /opt/local/include/eigen3 
)

IF(EIGEN_INCLUDE_DIRS)
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
  MESSAGE(STATUS "Found Eigen: ${EIGEN_INCLUDE_DIRS}")
ELSE(EIGEN_INCLUDE_DIRS)
  MESSAGE(FATAL_ERROR "Eigen not found")
ENDIF(EIGEN_INCLUDE_DIRS)

FIND_PATH(
  EIGEN_UN_INCLUDE_DIRS Eigen
  /usr/include/eigen3/unsupported
  /usr/local/include/eigen3/unsupported
  /opt/local/include/eigen3/unsupported
)

IF(EIGEN_UN_INCLUDE_DIRS)
  INCLUDE_DIRECTORIES(${EIGEN_UN_INCLUDE_DIRS})
  MESSAGE(STATUS "Found Unsupported Eigen: ${EIGEN_UN_INCLUDE_DIRS}")
ELSE(EIGEN_UN_INCLUDE_DIRS)
  MESSAGE(FATAL_ERROR "Unsupported Eigen not found")
ENDIF(EIGEN_UN_INCLUDE_DIRS)

FIND_PATH(
  UMFPACK_INCLUDE_DIRS umfpack.h
  /usr/include/suitesparse
  /usr/include/ufsparse
  /usr/local/include/suitesparse
  /usr/local/include/ufsparse
  /opt/local/include/suitesparse
  /opt/local/include/ufsparse
)

IF(UMFPACK_INCLUDE_DIRS)
  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
  ADD_DEFINITIONS(-DEIGEN_UMFPACK_SUPPORT)
  MESSAGE(STATUS "Found UMFPACK: ${UMFPACK_INCLUDE_DIRS}")
ELSE(UMFPACK_INCLUDE_DIRS)
  MESSAGE(FATAL_ERROR "UMFPACK not found")
ENDIF(UMFPACK_INCLUDE_DIRS)

