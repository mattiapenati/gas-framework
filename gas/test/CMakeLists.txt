ADD_EXECUTABLE(empty empty.cpp)
ADD_TEST(empty ${PROJECT_BINARY_DIR}/test/empty)

ADD_EXECUTABLE(gas_geometry_map_affine_interval geometry/map/affine_interval.cpp)
ADD_TEST(gas::geometry::map::affine<interval> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_interval)

ADD_EXECUTABLE(gas_geometry_map_affine_triangle geometry/map/affine_triangle.cpp)
ADD_TEST(gas::geometry::map::affine<triangle> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_triangle)

ADD_EXECUTABLE(gas_geometry_map_affine_triangle geometry/map/affine_triangle.cpp)
ADD_TEST(gas::geometry::map::affine<triangle> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_triangle)

ADD_EXECUTABLE(gas_numerical_quadrature_gauss_legendre_test numerical/quadrature/gauss_legendre.cpp)
ADD_TEST(gas::numerical::quadrature::gauss_legendre ${PROJECT_BINARY_DIR}/test/gas_numerical_quadrature_gauss_legendre_test)