FIND_PATH(EIGEN_INCLUDE_DIRS Eigen/Array
          /usr/include/eigen2
          /opt/local/include/eigen2 
         )

INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})