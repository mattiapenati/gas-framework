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
#include "Gas/Integration/Integration.h"
#include "Gas/Element/Element.h"

double u (double const & x, double const & y) {
	return 2.;
}
double up (double const & x, double const & y) {
	return 0.;
}

double f (double const & x, double const & y) {
	return u(x,y);
}

template <typename Return, typename Function>
Return eval (Function const & f, double const & x, double const & y) {
	return f(x,y);
}

#define PRINT_POL(a, b, c) \
	if (a == 0.) { \
		if ((b == c) && (c == 0.)) { \
			std::cout << "0."; \
		} else { \
			if (b != 0.) { \
				std::cout << b << "x"; \
				if (c != 0.) \
					std::cout << (c < 0 ? "" : "+") << c << "y"; \
			} else { \
				std::cout << c << "y"; \
			} \
		} \
	} else { \
		std::cout << a; \
		if (b != 0.) \
			std::cout << (b < 0 ? "" : "+") << b << "x"; \
		if (c != 0.) \
			std::cout << (c < 0 ? "" : "+") << c << "y"; \
	} \

struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

#define N 9

int main (int argc, char * argv[]) {

	std::cout.precision(3);

	/* matrice di stiffness, termine noto e soluzione */
	LinearAlgebra::Matrix<double, N, N> m(0.);
	LinearAlgebra::Vector<double, N> v(0.);
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
		++i;
		++itV;
	}
	std::cout << std::endl;

	/* integratore */
	typedef Integrator<Method::NewtonCotes<Geometry::Triangle, 1> > Integ;
	Integ integ;

	/* iterazione sulla faccia */
	typedef Triangulation::Finite_faces_iterator FaceIterator;
	FaceIterator itF = t.finite_faces_begin();
	int iF = 1;
	while (itF != t.finite_faces_end()) {
		Geometry::Triangle gt(*itF);
		Element::P1<Geometry::Triangle> el(gt);
		integ.domain(gt);

		/* indice dei nodi */
		LinearAlgebra::Vector<unsigned int, 3> index;
		for (unsigned int i = 0; i < 3; ++i) {
			index(i) = itF->vertex(i)->id();
		}

		/* operatori locali */
		for (unsigned int i = 0; i < 3; ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				m(index(i),index(j)) += integ.integrateMul<Integ::Transform, Integ::Transform>(el.grad(i), el.grad(j));
				m(index(i),index(j)) += integ.integrateMul<Integ::Transform, Integ::Transform>(el.base(i), el.base(j));
			}
		}

		/* termine noto locale */
		for (unsigned int i = 0; i < 3; ++i) {
			v(index(i)) += integ.integrateMul<Integ::Transform, Integ::Transform>(&f, el.base(i));
		}

		/* coordinate dei nodi */
		LinearAlgebra::Vector<double, 2> P0(itF->vertex(0)->point().x(), itF->vertex(0)->point().y());
		LinearAlgebra::Vector<double, 2> P1(itF->vertex(1)->point().x(), itF->vertex(1)->point().y());
		LinearAlgebra::Vector<double, 2> P2(itF->vertex(2)->point().x(), itF->vertex(2)->point().y());

		/* ricostruzione coefficienti delle funzioni a pezzi */
		LinearAlgebra::Matrix<double, 3, 3> M(1.);
		LinearAlgebra::Vector<double, 3> phi0, phi1, phi2, e;
		M(0,1) = P0(0); M(0,2) = P0(1);
		M(1,1) = P2(0); M(1,2) = P2(1);
		M(2,1) = P1(0); M(2,2) = P1(1);

		e(0) = eval<double>(el.base(0), P0(0), P0(1));
		e(1) = eval<double>(el.base(0), P2(0), P2(1));
		e(2) = eval<double>(el.base(0), P1(0), P1(1));
		LinearAlgebra::Solver::solver<LinearAlgebra::Solver::LU>::solve(M, phi0, e);
		e(0) = eval<double>(el.base(1), P0(0), P0(1));
		e(1) = eval<double>(el.base(1), P2(0), P2(1));
		e(2) = eval<double>(el.base(1), P1(0), P1(1));
		LinearAlgebra::Solver::solver<LinearAlgebra::Solver::LU>::solve(M, phi1, e);
		e(0) = eval<double>(el.base(2), P0(0), P0(1));
		e(1) = eval<double>(el.base(2), P2(0), P2(1));
		e(2) = eval<double>(el.base(2), P1(0), P1(1));
		LinearAlgebra::Solver::solver<LinearAlgebra::Solver::LU>::solve(M, phi2, e);
		
		/* ricostruzione gradienti */
		LinearAlgebra::Vector<double, 2> grad0 = eval<LinearAlgebra::Vector<double, 2> >(el.grad(0), 0., 0.);
		LinearAlgebra::Vector<double, 2> grad1 = eval<LinearAlgebra::Vector<double, 2> >(el.grad(1), 0., 0.);
		LinearAlgebra::Vector<double, 2> grad2 = eval<LinearAlgebra::Vector<double, 2> >(el.grad(2), 0., 0.);

		/* stampo informazioni */
		std::cout << "=== Faccia " << iF << " ===" << std::endl;
		std::cout << "P0(" << index(0) << "): " << P0(0) << " " << P0(1) << std::endl;
		std::cout << "P1(" << index(1) << "): " << P1(0) << " " << P1(1) << std::endl;
		std::cout << "P2(" << index(2) << "): " << P2(0) << " " << P2(1) << std::endl;
		std::cout << std::endl;
		std::cout << "phi0(x,y) = "; PRINT_POL(phi0(0), phi0(1), phi0(2)) std::cout << std::endl;
		std::cout << "phi1(x,y) = "; PRINT_POL(phi1(0), phi1(1), phi1(2)) std::cout << std::endl;
		std::cout << "phi2(x,y) = "; PRINT_POL(phi2(0), phi2(1), phi2(2)) std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "grad0(x,y) = [ " << grad0(0) << ", " << grad0(1) << " ]" << std::endl;
		std::cout << "grad1(x,y) = [ " << grad1(0) << ", " << grad1(1) << " ]" << std::endl;
		std::cout << "grad2(x,y) = [ " << grad2(0) << ", " << grad2(1) << " ]" << std::endl;
		std::cout << std::endl;

		++iF;
		++itF;
	}

	std::cout << "=== Matrice di stiffness ===" << std::endl;
	for (unsigned int i = 0; i < N; ++i) {
		for (unsigned int j = 0; j < N; ++j)
			std::cout << m(i,j) << "\t";
		std::cout << std::endl;
	}
	std::cout << std::endl;

	LinearAlgebra::Solver::solver<LinearAlgebra::Solver::LU>::solve(m, x, v);
	
	std::cout << "=== Soluzione ===" << std::endl;
	std::cout << "x = [ ";
	for (unsigned int i = 0; i < N; ++i)
		std::cout << x(i) << " ";
	std::cout << "]" << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Norma ===" << std::endl;
	double normL2 = 0.;
	double normH1 = 0.;
	{
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tL2 = el.normL2();
			double tH1 = el.normH1();

			normL2 += (tL2 * tL2);
			normH1 += (tH1 * tH1);
			
			++itF;
		}

		normL2 = std::sqrt(normL2);
		normH1 = std::sqrt(normH1);
	}
	std::cout << "Norma L2 = " << normL2 << std::endl;
	std::cout << "Norma H1 = " << normH1 << std::endl;
	std::cout << std::endl;
	
	std::cout << "=== Errore ===" << std::endl;
	double errL2 = 0.;
	double errH1 = 0.;
	{
		FaceIterator itF = t.finite_faces_begin();
		while (itF != t.finite_faces_end()) {
			Geometry::Triangle gt(*itF);
			Element::P1<Geometry::Triangle> el(gt);

			for (unsigned int i = 0; i < 3; ++i)
				el(i) = x(itF->vertex(i)->id());

			double tL2 = el.normL2(&u);
			double tH1 = el.normH1(&u, &up, &up);

			errL2 += (tL2 * tL2);
			errH1 += (tH1 * tH1);
			
			++itF;
		}

		errL2 = std::sqrt(errL2);
		errH1 = std::sqrt(errH1);
	}
	std::cout << "Errore L2 = " << errL2 << std::endl;
	std::cout << "Errore H1 = " << errH1 << std::endl;

	return 0;
}
