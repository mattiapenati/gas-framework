FIND_PATH(CGAL_DIR CGALConfig.cmake
          /usr/lib/CGAL
          /usr/local/lib/CGAL
          /opt/local/lib/CGAL
         )

IF(CGAL_DIR)
  SET(CGAL_FOUND TRUE)
  INCLUDE(${CGAL_DIR}/CGALConfig.cmake)
  INCLUDE(${CGAL_DIR}/UseCGAL.cmake)
ELSE(CGAL_DIR)
  MESSAGE(STATUS "CGAL not found.")
ENDIF(CGAL_DIR)
