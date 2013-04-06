
# Find header files
FIND_PATH(libconfig_INCLUDE_DIRS "libconfig.h++"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
)

# Find library
FIND_LIBRARY(libconfig_LIBRARIES
  NAMES "libconfig++.so"
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "/usr/local/lib"
)

# Check that we have found everything
SET(libconfig_FOUND FALSE)
IF(libconfig_INCLUDE_DIRS AND libconfig_LIBRARIES)
  SET(libconfig_FOUND TRUE)
ELSE(libconfig_INCLUDE_DIRS AND libconfig_LIBRARIES)
  MESSAGE(FATAL_ERROR "Libconfig not found - Please install Libconfig")
ENDIF(libconfig_INCLUDE_DIRS AND libconfig_LIBRARIES)
