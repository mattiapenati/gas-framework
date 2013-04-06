ADD_EXECUTABLE(img_seg main.cpp)
TARGET_LINK_LIBRARIES(img_seg X11 pthread umfpack amd blas ${DEPLIBS})
