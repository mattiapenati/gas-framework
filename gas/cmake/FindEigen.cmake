FIND_PATH(EIGEN_INCLUDE_DIRS Eigen/Array
          /usr/include/eigen2
          /usr/local/include/eigen2
          /opt/local/include/eigen2 
         )

IF(EIGEN_INCLUDE_DIRS)
  SET(EIGEN_FOUND TRUE)
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
  MESSAGE(STATUS "Found Eigen: ${EIGEN_INCLUDE_DIRS}")
ENDIF(EIGEN_INCLUDE_DIRS)

FIND_PATH(UMFPACK_INCLUDE_DIRS umfpack.h
          /usr/include/suitesparse
          /usr/include/ufsparse
          /usr/local/include/suitesparse
          /usr/local/include/ufsparse
          /opt/local/include/suitesparse
          /opt/local/include/ufsparse
         )
IF(UMFPACK_INCLUDE_DIRS)
  SET(UMFPACK_FOUND TRUE)
  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
  ADD_DEFINITIONS(-DEIGEN_UMFPACK_SUPPORT)
  MESSAGE(STATUS "Found UmfPack: ${UMFPACK_INCLUDE_DIRS}")
ENDIF(UMFPACK_INCLUDE_DIRS)
          