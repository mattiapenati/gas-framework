#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_vertex_base_with_id_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <iostream>
#include <cmath>

#include "Gas/Geometry/Geometry.h"
#include "Gas/LinearAlgebra/LinearAlgebra.h"
#include "Gas/Element/Element.h"

double u (double const & x, double const & y) {
	return 2. + 3. * x + 4. * y;
}
double ux (double const & x, double const & y) {
	return 3.;
}
double uy (double const & x, double const & y) {
	return 4.;
}

double f (double const & x, double const & y) {
	return u(x,y);
}

template <typename Return, typename Function>
Return eval (Function const & f, double const & x, double const & y) {
	return f(x,y);
}

struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

#define N 9

int main (int argc, char * argv[]) {

	std::cout.precision(3);

	/* soluzione */
	LinearAlgebra::Vector<double, N> x(0.);

	/* creating a triangulation */
	typedef CGAL::Triangulation_vertex_base_with_id_2<K> Vb;
	typedef CGAL::Delaunay_mesh_face_base_2<K, CGAL::Constrained_triangulation_face_base_2<K> > Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> Triangulation;
	typedef CGAL::Delaunay_mesh_size_criteria_2<Triangulation> Criteria;
	Triangulation t;

	{
		/* inserimento punti */
		typedef Triangulation::Point Point;
		Point a(+1.,+1.);
		Point b(-1.,+1.);
		Point c(-1.,-1.);
		Point d(+1.,-1.);
	
		Triangulation::Vertex_handle va = t.insert(a);
		Triangulation::Vertex_handle vb = t.insert(b);
		Triangulation::Vertex_handle vc = t.insert(c);
		Triangulation::Vertex_handle vd = t.insert(d);
	
		/* inserimento lati */
		t.insert_constraint(va, vb);
		t.insert_constraint(vb, vc);
		t.insert_constraint(vc, vd);
		t.insert_constraint(vd, va);
	
		/* raffinamento */
		CGAL::refine_Delaunay_mesh_2(t, Criteria(0.125, 1.5));
	}

	/* numerazione dei punti */
	std::cout << "=== Numerazione nodi ===" << std::endl;
	typedef Triangulation::Finite_vertices_iterator VerticesIterator;
	int i = 0;
	VerticesIterator itV = t.finite_vertices_begin();
	while (itV != t.finite_vertices_end()) {
		std::cout << "P"  << i << " : " << itV->point().x() << " " << itV->point().y() << std::endl;
		itV->id() = i;
		x(i) = f(itV->point().x(), itV->point().y());
		++i;
		++itV;
	}
	std::cout << std::endl;
	
	std::cout << "=== Soluzione ===" << std::endl;
	std::cout << "x = [ ";
	for (unsigned int i = 0; i < N; ++i)
		std::cout << x(i) << " ";
	std::cout << "]" << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Norma L2 ===" << std::endl;
	double normL2 = 0.;
	{
		typedef Triangulation::Finite_faces_iterator FaceIterator;
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tL2 = el.normL2(); tL2 *= tL2;
			std::cout << tL2 << "\t";
			normL2 += tL2;
			
			++itF;
		}

		normL2 = std::sqrt(normL2);
	}
	std::cout << std::endl;
	std::cout << "Norma L2 = " << normL2 << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Norma H1 ===" << std::endl;
	double normH1 = 0.;
	{
		typedef Triangulation::Finite_faces_iterator FaceIterator;
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tH1 = el.normH1(); tH1 *= tH1;
			std::cout << tH1 << "\t";
			normH1 += tH1;
			
			++itF;
		}

		normH1 = std::sqrt(normH1);
	}
	std::cout << std::endl;
	std::cout << "Norma H1 = " << normH1 << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Errore L2 ===" << std::endl;
	double errL2 = 0.;
	{
		typedef Triangulation::Finite_faces_iterator FaceIterator;
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tL2 = el.normL2(&u); tL2 *= tL2;
			std::cout << tL2 << "\t";
			errL2 += tL2;
			
			++itF;
		}

		errL2 = std::sqrt(errL2);
	}
	std::cout << std::endl;
	std::cout << "Errore L2 = " << errL2 << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Errore H1 ===" << std::endl;
	double errH1 = 0.;
	{
		typedef Triangulation::Finite_faces_iterator FaceIterator;
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tH1 = el.normH1(&u, &ux, &uy); tH1 *= tH1;
			std::cout << tH1 << "\t";
			errH1 += tH1;
			
			++itF;
		}

		errH1 = std::sqrt(errH1);
	}
	std::cout << std::endl;
	std::cout << "Errore H1 = " << errH1 << std::endl;
	std::cout << std::endl;

	return 0;
}
