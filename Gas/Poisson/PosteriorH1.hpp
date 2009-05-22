#ifndef GAS_POSTERIORH1_HPP
#define GAS_POSTERIORH1_HPP

#include "../Integration/Integration.h"
#include "../LinearAlgebra/LinearAlgebra.h"
#include "../Element/Element.h"

#ifndef MAX
#define MAX(A, B) ((A > B) ? A : B)
#endif

template < typename Triangulation >
class PosteriorH1 {

	private:
		typedef PosteriorH1<Triangulation> self;
		
		// Integratore
		typedef Integrator < Method::NewtonCotes<Geometry::Triangle, 2> >  Integrator2;
		
		// Algebra lineare
		typedef LinearAlgebra::Vector<double, 2> vector2;
		typedef LinearAlgebra::Matrix<double, 2, 2> matrix22;

		// Elementi
		typedef Element::P1<Geometry::Triangle> ElementP1;

	public:
		PosteriorH1 (Triangulation const &);
		
		template < typename FacePointer, typename Function >
		double residue (FacePointer, Function const &);
	
	private:
		// Base P1
		static double P1_phi0 (double const & x, double const & y) { return 1 - x -y; }
		static double P1_phi1 (double const & x, double const & y) { return x; }
		static double P1_phi2 (double const & x, double const & y) { return y; }
		
		// Base Clement
		static double C_phi0 (double const & x , double const & y) { return 1.; }
		static double C_phi1 (double const & x , double const & y) { return x; }
		static double C_phi2 (double const & x , double const & y) { return y; }
	
	private:
		// Costruzione del gradiente
		vector2 grad(vector2 const &, double const &, vector2 const &, double const &, vector2 const &, double const &);
		// Costruzione interpolante di Clement in un nodo (restituisce il gradiente)
		template < typename VertexPointer >
		vector2 grad_Clement(VertexPointer);
		
	private:
		Triangulation const & t_;
		
		double normH1;

};

/*
 * Costruttore
 */
template <typename Triangulation >
PosteriorH1<Triangulation>::PosteriorH1 (Triangulation const & t): t_(t), normH1(0.) {
	// Iteratore su tutte le facce della triangolazione
	typename Triangulation::Finite_faces_iterator itF = t_.finite_faces_begin();
	while (itF != t_.finite_faces_end()) {
		// Soluzione
		LinearAlgebra::Vector<double, 3> u;
		for (unsigned int i = 0; i < 3; ++i)
			u(i) = itF->vertex(i)->info().value();

		// Elemento
		Geometry::Triangle g(*itF);
		ElementP1 e(g, u);

		double t = e.normH1();
		normH1 += (t*t);

		// Avanzamento dell'iteratore
		++itF;
	}
	normH1 = std::sqrt(normH1);
}

/*
 * Calcolo del residuo
 */
template <typename Triangulation>
template <typename FacePointer, typename Function >
double PosteriorH1<Triangulation>::residue (FacePointer f, Function const & forzante) {
	
	double res = 0.;
	
	// Coordinate dei nodi
	LinearAlgebra::Vector<double, 2> p0(f->vertex(0)->point().x(), f->vertex(0)->point().y());
	LinearAlgebra::Vector<double, 2> p1(f->vertex(1)->point().x(), f->vertex(1)->point().y());
	LinearAlgebra::Vector<double, 2> p2(f->vertex(2)->point().x(), f->vertex(2)->point().y());
	
	// Distanze
	double d01 = norm(p0-p1);
	double d02 = norm(p0-p2);
	double d12 = norm(p1-p2);
	double h = MAX(d01, MAX(d02, d12));
	
	// Soluzione
	LinearAlgebra::Vector<double, 3> u;
	u(0) = f->vertex(0)->info().value();
	u(1) = f->vertex(1)->info().value();
	u(2) = f->vertex(2)->info().value();
	
	// Elemento
	Geometry::Triangle g(*f);
	Element::P1<Geometry::Triangle> e(g, u);
	
	// Gradiente locale
	LinearAlgebra::Vector<double, 2> grad_loc = e.grad(0.,0.);
	
	// Gradiente di Clement nei tre vertici
	LinearAlgebra::Vector<double, 2> gradC0 = grad_Clement(f->vertex(0));
	LinearAlgebra::Vector<double, 2> gradC1 = grad_Clement(f->vertex(1));
	LinearAlgebra::Vector<double, 2> gradC2 = grad_Clement(f->vertex(2));
	
	// Laplaciano
	LinearAlgebra::Vector<double, 2> lap0 = grad(p0, gradC0(0), p1, gradC1(0), p2, gradC2(0));
	LinearAlgebra::Vector<double, 2> lap1 = grad(p0, gradC0(1), p1, gradC1(1), p2, gradC2(1));
	double laplaciano = lap0(0) + lap1(1);
	
	// Integratore
	Integrator2 integr(g);
	
	// Residuo sulla faccia
	double res_faccia = h * std::sqrt(
			integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(forzante , forzante) + 
			2 * laplaciano * integr.integrate<Integrator2::Transform>(forzante) +
			laplaciano * laplaciano * g.area()
		);
	
	// Residuo lato
	double res_lato = 0.;
	// Giro sui lati
	for (unsigned int i = 0 ; i < 3 ; ++i) {
		
		// Controllo se faccia finita
		if (t_.is_infinite(f->neighbor(i)))
			continue;
		
		// Coordinate dei nodi del vicino
		LinearAlgebra::Vector<double, 2> pn0(f->neighbor(i)->vertex(0)->point().x(), f->neighbor(i)->vertex(0)->point().y());
		LinearAlgebra::Vector<double, 2> pn1(f->neighbor(i)->vertex(1)->point().x(), f->neighbor(i)->vertex(1)->point().y());
		LinearAlgebra::Vector<double, 2> pn2(f->neighbor(i)->vertex(2)->point().x(), f->neighbor(i)->vertex(2)->point().y());
	
		// Soluzione del vicino
		LinearAlgebra::Vector<double, 3> un;
		un(0) = f->neighbor(i)->vertex(0)->info().value();
		un(1) = f->neighbor(i)->vertex(1)->info().value();
		un(2) = f->neighbor(i)->vertex(2)->info().value();
		
		// Gradiente del vicino
		LinearAlgebra::Vector<double, 2> grad_n = grad(pn0, un(0), pn1, un(1), pn2, un(2));
		
		// Normale
		LinearAlgebra::Vector<double, 2> n;
		double t;
		switch(i) {
			case 0:
				n(0) = p1(1) - p2(1);
				n(1) = p2(0) - p1(0);
				t = dot(grad_loc - grad_n, n);
				res_lato += t * t / d12;
				break;
			case 1:
				n(0) = p2(1) - p0(1);
				n(1) = p0(0) - p2(0);
				t = dot(grad_loc - grad_n, n);
				res_lato += t * t / d02;
				break;
			case 2:
				n(0) = p0(1) - p1(1);
				n(1) = p1(0) - p0(0);
				t = dot(grad_loc - grad_n, n);
				res_lato += t * t / d01;
				break;
		}
	}
	
	res_lato = std::sqrt(h * res_lato) / 2;
	
	res = res_faccia + res_lato;
	
	return res / normH1;
}

/*
 * Calcolo del gradiente di una funzione P1, partendo dai tre valori nei nodi
 */
template < typename Triangulation >
LinearAlgebra::Vector<double, 2> PosteriorH1<Triangulation>::grad(
		vector2 const & p0, double const & u0,
		vector2 const & p1, double const & u1,
		vector2 const & p2, double const & u2)
{
	vector2 gradr;
	
	// Distanze
	double d01 = norm(p0-p1);
	double d02 = norm(p0-p2);
	double d12 = norm(p1-p2);
	
	// Coefficienti del gradiente
	double a01 = (u1 - u0) / (d01 * d01);
	double a02 = (u2 - u0) / (d02 * d02);
	
	gradr = (p1-p0) * a01 + (p2-p0)*a02;
	
	return gradr;
}
/*
 * Costruzione interpolante di Clement in un nodo (restituisce il gradiente)
 */
template < typename Triangulation>
template < typename VertexPointer >
LinearAlgebra::Vector<double, 2> PosteriorH1<Triangulation>::grad_Clement (VertexPointer v) {
	
	// Coordinate del punto
	double x = v->point().x();
	double y = v->point().y();
		
	// Matrice di stiffness
	LinearAlgebra::Matrix<double,3,3> M(0.);
	// Vettore termine noto
	LinearAlgebra::Vector<double,3> bx(0.) , by(0.);
	// Vettore soluzione
	LinearAlgebra::Vector<double,3> ax(0.) , ay(0.);
	
	// Integratore
	Integrator2 integr;
	
	// Gradiente di ritorno
	vector2 gradC;
	
	// Ciclo sulle facce incidenti al vertice
	typename Triangulation::Face_circulator cF, end;
	cF = t_.incident_faces(v);
	end = cF;
	do {
		// Matrice di proiezione locale
		if (!t_.is_infinite(cF)) {
			Integrator2::Geometry g(*cF);
			integr.domain(g);
			
			// Matrice
			M(0,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi0, self::C_phi0);
			M(0,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi0, self::C_phi1);
			M(0,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi0, self::C_phi2);
			M(1,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi1, self::C_phi0);
			M(1,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi1, self::C_phi1);
			M(1,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi1, self::C_phi2);
			M(2,0) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi2, self::C_phi0);
			M(2,1) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi2, self::C_phi1);
			M(2,2) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(self::C_phi2, self::C_phi2);
			
			// Coordinate dei nodi
			LinearAlgebra::Vector<double, 2> p0(cF->vertex(0)->point().x(), cF->vertex(0)->point().y());
			LinearAlgebra::Vector<double, 2> p1(cF->vertex(1)->point().x(), cF->vertex(1)->point().y());
			LinearAlgebra::Vector<double, 2> p2(cF->vertex(2)->point().x(), cF->vertex(2)->point().y());
			
			// Soluzione
			LinearAlgebra::Vector<double, 3> u;
			u(0) = cF->vertex(0)->info().value();
			u(1) = cF->vertex(1)->info().value();
			u(2) = cF->vertex(2)->info().value();
			
			// Gradiente locale
			LinearAlgebra::Vector<double, 2> grad_loc = grad(p0, u(0), p1, u(1), p2, u(2));
			
			// Termini noti
			bx(0) += grad_loc(0) * integr.integrate<Integrator2::Transform>(self::C_phi0);
			bx(1) += grad_loc(0) * integr.integrate<Integrator2::Transform>(self::C_phi1);
			bx(2) += grad_loc(0) * integr.integrate<Integrator2::Transform>(self::C_phi2);
			
			by(0) += grad_loc(1) * integr.integrate<Integrator2::Transform>(self::C_phi0);
			by(1) += grad_loc(1) * integr.integrate<Integrator2::Transform>(self::C_phi1);
			by(2) += grad_loc(1) * integr.integrate<Integrator2::Transform>(self::C_phi2);
		}
		++cF;
	} while(cF != end);
	
	// Soluzione del sistema
	LinearAlgebra::Solver::LU::solve( M , ax , bx );
	LinearAlgebra::Solver::LU::solve( M , ay , by );
	
	// Gradiente di Clement
	gradC(0) = ax(0) + ax(1) * x + ax(2) * y;
	gradC(1) = ay(0) + ay(1) * x + ay(2) * y;
	
	return gradC;
}

#endif // GAS_POSTERIORH1_HPP