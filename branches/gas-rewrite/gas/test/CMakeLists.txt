ADD_EXECUTABLE(empty empty.cpp)
ADD_TEST(empty ${PROJECT_BINARY_DIR}/test/empty)

ADD_EXECUTABLE(gas_functional_base_P1_interval functional/base/P1_interval.cpp)
ADD_TEST(gas::functional::base::P1<interval> ${PROJECT_BINARY_DIR}/test/gas_functional_base_P1_interval)

ADD_EXECUTABLE(gas_functional_base_P1_triangle functional/base/P1_triangle.cpp)
ADD_TEST(gas::functional::base::P1<triangle> ${PROJECT_BINARY_DIR}/test/gas_functional_base_P1_triangle)

ADD_EXECUTABLE(gas_functional_base_P1_tetra functional/base/P1_tetra.cpp)
ADD_TEST(gas::functional::base::P1<tetra> ${PROJECT_BINARY_DIR}/test/gas_functional_base_P1_tetra)

ADD_EXECUTABLE(gas_geometry_unit geometry/unit.cpp)
ADD_TEST(gas::geometry::unit ${PROJECT_BINARY_DIR}/test/gas_geometry_unit)

ADD_EXECUTABLE(gas_geometry_map_affine_interval geometry/map/affine_interval.cpp)
ADD_TEST(gas::geometry::map::affine<interval> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_interval)

ADD_EXECUTABLE(gas_geometry_map_affine_triangle geometry/map/affine_triangle.cpp)
ADD_TEST(gas::geometry::map::affine<triangle> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_triangle)

ADD_EXECUTABLE(gas_geometry_map_affine_tetra geometry/map/affine_tetra.cpp)
ADD_TEST(gas::geometry::map::affine<tetra> ${PROJECT_BINARY_DIR}/test/gas_geometry_map_affine_tetra)

#ADD_EXECUTABLE(gas_numerical_quadrature_gauss_legendre_test numerical/quadrature/gauss_legendre.cpp)
#ADD_TEST(gas::numerical::quadrature::gauss_legendre ${PROJECT_BINARY_DIR}/test/gas_numerical_quadrature_gauss_legendre_test)

#ADD_EXECUTABLE(gas_numerical_quadrature_newtoncotes_test numerical/quadrature/newton_cotes.cpp)
#ADD_TEST(gas::numerical::quadrature::newton_cotes ${PROJECT_BINARY_DIR}/test/gas_numerical_quadrature_newtoncotes_test)