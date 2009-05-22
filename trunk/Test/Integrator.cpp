#include "Gas/Geometry/Geometry.h"
#include "Gas/Integration/Integration.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include <iostream>

double f (double const & x, double const & y) {
	return x;
}
double g (double const & x, double const & y) {
	return x;
}
double fg (double const & x, double const & y) {
	return x*x;
}
double uno (double const & x, double const & y) {
	return 1.;
}

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

	/* integratore */
	typedef Integrator<Method::NewtonCotes<Geometry::Triangle, 1> > Integ;
	Integ integ;
	double s = 0.;

	/* iterazione sulla faccia */
	typedef Triangulation::Finite_faces_iterator FaceIterator;
	FaceIterator itF = t.finite_faces_begin();
	while (itF != t.finite_faces_end()) {
		Geometry::Triangle gt(*itF);
		integ.domain(gt);

		double mi = integ.integrateMul<Integ::Transform, Integ::Transform>(&f, &g);
		double i  = integ.integrate<Integ::Transform>(&fg);
		
		s += mi;

		std::cout << mi << " = " << i << std::endl;

		++itF;
	}
	std::cout << s << std::endl;

	/* geometria di default */
	Geometry::Triangle gt;
	
	integ.domain(gt);

	std::cout << gt.area() << " = " << integ.integrate<Integ::Transform>(&uno) << std::endl;

	return 0;
}
