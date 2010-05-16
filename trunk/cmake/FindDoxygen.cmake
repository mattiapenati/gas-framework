FIND_PROGRAM(
  DOXYGEN_COMMAND
  doxygen
  PATHS /usr/bin /usr/local/bin /opt/local/bin
)

IF(DOXYGEN_COMMAND)
  SET(DOXYGEN_FOUND TRUE)
  MESSAGE(STATUS "Found Doxygen: ${DOXYGEN_COMMAND}")
  MARK_AS_ADVANCED(DOXYGEN_FOUND)
ENDIF(DOXYGEN_COMMAND)
  