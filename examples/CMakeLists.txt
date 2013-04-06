if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


#ADD_EXECUTABLE(empty empty.cpp)

#ADD_SUBDIRECTORY(poisson)

ADD_SUBDIRECTORY(image_segmentation)
