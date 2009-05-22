add_executable(Triangle Triangle.cpp)
add_executable(Integrator Integrator.cpp)
add_executable(Element Element.cpp)
add_executable(Norm Norm.cpp)

target_link_libraries(Triangle CGAL)
target_link_libraries(Integrator CGAL)
target_link_libraries(Element CGAL)
target_link_libraries(Norm CGAL)
