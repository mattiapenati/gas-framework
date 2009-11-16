IF(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
  SET(CGAL_FOUND TRUE)
ELSE(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
  FIND_PATH(CGAL_DIR CGALConfig.cmake
            /usr/lib/CGAL
            /usr/local/lib/CGAL
            /opt/local/lib/CGAL
           )
  IF(CGAL_DIR)
    SET(CGAL_FOUND TRUE)
    INCLUDE(${CGAL_DIR}/CGALConfig.cmake)
    INCLUDE(${CGAL_DIR}/UseCGAL.cmake)
    INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
    LINK_DIRECTORIES(${CGAL_LIBRARIES_DIR})
    MESSAGE(STATUS "Found CGAL: ${CGAL_INCLUDE_DIRS}, ${CGAL_LIBRARIES_DIR}")
  ELSE(CGAL_DIR)
    MESSAGE(STATUS "CGAL not found.")
  ENDIF(CGAL_DIR)
ENDIF(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)