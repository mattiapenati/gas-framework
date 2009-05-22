#include "Gas/Geometry/Geometry.h"
#include "Gas/LinearAlgebra/LinearAlgebra.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include <iostream>

struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

int main (int argc, char * argv[]) {

	/* creating a triangulation */
	typedef CGAL::Triangulation_2<K> Triangulation;
	Triangulation t;

	/* inserimento punti */
	typedef Triangulation::Point Point;
	Point a(+1.,+1.);
	Point b(-1.,+1.);
	Point c(-1.,-1.);
	Point d(+1.,-1.);

	t.insert(a);
	t.insert(b);
	t.insert(c);
	t.insert(d);

	std::cout << "TRIANGOLAZIONE" << std::endl << std::endl;

	/* iterazione sulla faccia */
	typedef Triangulation::Finite_faces_iterator FaceIterator;
	FaceIterator itF = t.finite_faces_begin();
	while (itF != t.finite_faces_end()) {
		Geometry::Triangle g(*itF);

		LinearAlgebra::Vector<double, 2> p0(g.x(0.,0.), g.y(0.,0.));
		LinearAlgebra::Vector<double, 2> p1(g.x(1.,0.), g.y(1.,0.));
		LinearAlgebra::Vector<double, 2> p2(g.x(0.,1.), g.y(0.,1.));

		std::cout << "=== Faccia ===" << std::endl;
		std::cout << " P0: " << p0(0) << ", " << p0(1) << std::endl;
		std::cout << " P1: " << p1(0) << ", " << p1(1) << std::endl;
		std::cout << " P2: " << p2(0) << ", " << p2(1) << std::endl;

		std::cout << "=== Inversion ===" << std::endl;
		std::cout << " P0: " << g.invx(p0(0), p0(1)) << ", " << g.invy(p0(0), p0(1)) << std::endl;
		std::cout << " P1: " << g.invx(p1(0), p1(1)) << ", " << g.invy(p1(0), p1(1)) << std::endl;
		std::cout << " P2: " << g.invx(p2(0), p2(1)) << ", " << g.invy(p2(0), p2(1)) << std::endl;

		std::cout << "=== Dati ===" << std::endl;
		std::cout << " Area: " << g.area() << std::endl;
		std::cout << " Det: " << g.det(0.,0.) << std::endl;

		std::cout << std::endl;

		++itF;
	}

	/* geometria di default */
	std::cout << "DEFAULT" << std::endl << std::endl;

	Geometry::Triangle g;

	LinearAlgebra::Vector<double, 2> p0(g.x(0.,0.), g.y(0.,0.));
	LinearAlgebra::Vector<double, 2> p1(g.x(1.,0.), g.y(1.,0.));
	LinearAlgebra::Vector<double, 2> p2(g.x(0.,1.), g.y(0.,1.));

	std::cout << "=== Faccia ===" << std::endl;
	std::cout << " P0: " << p0(0) << ", " << p0(1) << std::endl;
	std::cout << " P1: " << p1(0) << ", " << p1(1) << std::endl;
	std::cout << " P2: " << p2(0) << ", " << p2(1) << std::endl;

	std::cout << "=== Inversion ===" << std::endl;
	std::cout << " P0: " << g.invx(p0(0), p0(1)) << ", " << g.invy(p0(0), p0(1)) << std::endl;
	std::cout << " P1: " << g.invx(p1(0), p1(1)) << ", " << g.invy(p1(0), p1(1)) << std::endl;
	std::cout << " P2: " << g.invx(p2(0), p2(1)) << ", " << g.invy(p2(0), p2(1)) << std::endl;

	std::cout << "=== Dati ===" << std::endl;
	std::cout << " Area: " << g.area() << std::endl;
	std::cout << " Det: " << g.det(0.,0.) << std::endl;

	return 0;
}
