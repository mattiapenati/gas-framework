FIND_PATH(CGAL_DIR CGALConfig.cmake
          /usr/lib/CGAL
          /usr/lib64/CGAL
          /usr/local/lib/CGAL
          /usr/local/lib64/CGAL
          /opt/local/lib/CGAL
         )

SET(CGAL_FOUND FALSE)
IF(CGAL_DIR)
  SET(CGAL_FOUND TRUE)
  INCLUDE(${CGAL_DIR}/CGALConfig.cmake)
  INCLUDE(${CGAL_DIR}/UseCGAL.cmake)
ELSE(CGAL_DIR)
  MESSAGE(FATAL ERROR "CGAL not found - Please install CGAL")
ENDIF(CGAL_DIR)
