#ifndef POISSON_HPP
#define POISSON_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "Info.hpp"
#include "PosteriorH1.hpp"

#include "../SparseLib++/coord_double.h"
#include "../SparseLib++/comprow_double.h"
#include "../SparseLib++/diagpre_double.h"
#include "../Iml/cg.h"

#include "../Integration/Integration.h"
#include "../LinearAlgebra/LinearAlgebra.h"
#include "../Element/Element.h"

#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#ifndef MAX
#define MAX(A, B) ((A > B) ? A : B)
#endif

#ifndef MIN
#define MIN(A, B) ((A < B) ? A : B)
#endif

class Poisson {
	private:
		// Kernel
		struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
		
		// Triangolazione
		typedef CGAL::Triangulation_vertex_base_with_info_2<PointInfo, K> Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K, 
			CGAL::Constrained_triangulation_face_base_2<K, 
				CGAL::Triangulation_face_base_with_info_2<FaceInfo , K>
			>
		> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	
	public:
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
	
	private:
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
		
		typedef double (*func)(double const &, double const &);
	public:
		typedef CDT::Point Point;
		typedef std::list<Point> PointList;

	public:
		Poisson(PointList const &, func, double const & criteria = 0.);
		~Poisson();
		template < typename Stimator > unsigned int refine(double const &);
		
		/* errori */
		template <typename Function> double errL2(Function const &);
		
		template <typename Function1, typename Function2, typename Function3>
		double errH1(Function1 const &, Function2 const &, Function3 const &);
		
		/* norme */
		double normL2();
		double normH1();
		
		void saveToSVG(char const *, bool const);
	private:
		// LinearAlgebra
		typedef CompRow_Mat_double Matrix;
		typedef MV_Vector_double Vector;
		
		// Integrator
		typedef Geometry::Triangle Triangle;
		typedef Integrator<Method::NewtonCotes<Triangle, 2> > Integrator2;

		// Element
		typedef Element::P1<Triangle> ElementP1;

		CDT cdt_;
		
		unsigned int n_point;
		
		func forzante;
		
		double area_minima;
	private:
		// Base P1
		static double P1_phi0 (double const & x, double const & y) { return 1 - x -y; }
		static double P1_phi1 (double const & x, double const & y) { return x; }
		static double P1_phi2 (double const & x, double const & y) { return y; }
		
		// Generazione griglia
		void insertNodes(PointList const &, double const &);
		void enumNodes();
		void setBConditions();
		
		// Costruzione matrice
		void makeMatrixTermineNoto(Matrix &, Vector &);
		void printPattern(char const *, Matrix const &);
		
		// Risoluzione del sistema lineare
		void solveSystem(Matrix const &, Vector &, Vector const &);
		
		// Salvataggio soluzione
		void saveSolution(Vector const &);
		
		// Calcolo dei residui
		void calcRes();
		
};

/*
 * Costruttore di default
 */
Poisson::Poisson(PointList const &boundary, func f, double const & criteria): forzante(f) {
	Matrix A;
	Vector x;
	Vector F;

	// Costruzione della griglia
	insertNodes(boundary, criteria);
	setBConditions();
	enumNodes();

	// Costruzione matrice termine noto
	makeMatrixTermineNoto(A, F);

	// Soluzione del sistema lineare
	solveSystem(A, x, F);

	// Salvataggio della soluzione nella griglia
	saveSolution(x);
}

/*
 * Distruttore
 */
Poisson::~Poisson() {
}

/*
 * Salva in un file SVG
 */
void Poisson::saveToSVG(char const *filename, bool const triangulation=true) {
	int H = 500.;
	int W = 500.;
	std::ofstream f;
	f.open(filename);
	f  << "<?xml version=\"1.0\" standalone=\"yes\"?>\n" \
		<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n" \
		<< "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	// Calcolo altezza e larghezza
	double xmin, xmax, ymin, ymax, umin, umax, rmin, rmax;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	xmin = xmax = itV->point().x();
	ymin = ymax = itV->point().y();
	umin = umax = itV->info().value();
	++itV;
	while(itV != cdt_.finite_vertices_end()) {
		double x = itV->point().x();
		double y = itV->point().y();
		double u = itV->info().value();
		xmin = MIN(xmin, x);
		xmax = MAX(xmax, x);
		ymin = MIN(ymin, y);
		ymax = MAX(ymax, y);
		umin = MIN(umin, u);
		umax = MAX(umax, u);
		++itV;
	}
	W = H * (xmax -xmin) / (ymax - ymin);
	// Inizio SVG
	f << "<svg width=\""<<W<<"px\" height=\""<<H<<"px\" " \
	  << "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	// Triangolazione
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	rmin = itF->info().res();
	rmax = itF->info().res();
	++itF;
	while(itF != cdt_.finite_faces_end()) {
		double res = itF->info().res();
		rmax = MAX(rmax, res);
		rmin = MIN(rmin, res);
		++itF;
	}
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Calcolo dei punti
		double x1 = itF->vertex(0)->point().x(); double y1 = itF->vertex(0)->point().y();
		double x2 = itF->vertex(1)->point().x(); double y2 = itF->vertex(1)->point().y();
		double x3 = itF->vertex(2)->point().x(); double y3 = itF->vertex(2)->point().y();
		x1 = (x1 - xmin) / (xmax - xmin); y1 = 1 - (y1 - ymin) / (ymax - ymin);
		x2 = (x2 - xmin) / (xmax - xmin); y2 = 1 - (y2 - ymin) / (ymax - ymin);
		x3 = (x3 - xmin) / (xmax - xmin); y3 = 1 - (y3 - ymin) / (ymax - ymin);
		x1 = x1 * W; y1 = y1 * H;
		x2 = x2 * W; y2 = y2 * H;
		x3 = x3 * W; y3 = y3 * H;
		
		double u = 0.;
		u += itF->vertex(0)->info().value();
		u += itF->vertex(1)->info().value();
		u += itF->vertex(2)->info().value();
		u /= 3;
		double res = itF->info().res();
		
		// Calcolo dei colori
		int color = 255 * (u - umin ) / (umax - umin);
		// color = 255 * (res - rmin ) / (rmax - rmin);
		int r = color;
		int g = color;
		int b = 0;
		
		// Stampa a file
		if (triangulation)
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"white\" width=\"1\">\n";
		else
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"rgb("<<r<<", "<<g<<", "<<b<<")\" width=\"1\">\n";
		f << "\t<polygon points=\""<<static_cast<int>(x1)<<", "<<static_cast<int>(y1)<<" "\
			<< static_cast<int>(x2)<<", "<<static_cast<int>(y2)<<" "\
			<< static_cast<int>(x3)<<", "<<static_cast<int>(y3)<<"\" />\n";
		f << "</g>\n";
	}
	f  << "</svg>\n";
}

/*
 * Pattern della matrice
 */
void Poisson::printPattern(char const * name, Matrix const & A) {
	std::fstream pbm;
	pbm.open(name, std::ios_base::out);
	
	unsigned int cols = A.dim(0);
	unsigned int rows = A.dim(1);
	
	pbm<<"P1"<<std::endl<<cols<<" "<<rows<<std::endl;
	
	for (unsigned int i = 0 ; i < rows; ++i) {
		for (unsigned int j = 0 ; j < cols; ++j) {
			if (A(i, j) == 0)
				pbm<<"0 ";
			else
				pbm<<"1 ";
			pbm<<std::endl;
		}
	}
}

/*
 * Inserimento dei nodi
 */
void Poisson::insertNodes(PointList const &boundary_, double const & criteria) {
	
	// Criterio mesher
	double minimo = 0.;
	double d;
	double x0, y0, x1, y1;
	
	// Vertex handle list
	typedef std::list<CDT::Vertex_handle> VertexList;
	VertexList Pb_;
	
	// Inserimento punti
	PointList::const_iterator itP;
	
	for (itP = boundary_.begin(); itP != boundary_.end(); ++itP)
		Pb_.push_back(cdt_.insert(*itP));
	
	// Inserimento lati constrained
	VertexList::const_iterator itV1 = Pb_.begin();
	VertexList::const_iterator itV2 = Pb_.begin();
	++itV2;
	
	LinearAlgebra::Vector<double, 2> p1((*itV1)->point().x(), (*itV1)->point().y());
	LinearAlgebra::Vector<double, 2> p2((*itV2)->point().x(), (*itV2)->point().y());
	
	minimo = norm(p1-p2);
	
	while(itV2 != Pb_.end()) {
		cdt_.insert_constraint(*itV1, *itV2);
		
		// lunghezza del lato
		p1(0) = (*itV1)->point().x(); p1(1) = (*itV1)->point().y();
		p2(0) = (*itV2)->point().x(); p2(1) = (*itV2)->point().y();
	
		d = norm(p1-p2);
		d = MIN(d, minimo);
		
		// avanzamento
		++itV1;
		++itV2;
	}
	
	// Ultimo lato
	cdt_.insert_constraint(Pb_.front(), Pb_.back());
	
	p1(0) = Pb_.front()->point().x(); p1(1) = Pb_.front()->point().y();
	p2(0) = Pb_.back()->point().x(); p2(1) = Pb_.back()->point().y();

	d = norm(p1-p2);
	d = MIN(d, minimo);
	
	if (criteria) minimo = criteria;
	
	// Impostazione dei nodi della lista al bordo (da non eliminare)
	CDT::Vertex_circulator cV = cdt_.incident_vertices (cdt_.infinite_vertex());
	CDT::Vertex_circulator begin = cV;
	do {
		cV->info().toBoundary();
		++cV;
	} while(cV != begin);
	
	// Raffinamento
	CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, minimo));
	
	// Calcolo dell'area minima
	{
		area_minima = 0.;
		CDT::Finite_faces_iterator itF;
		for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
			Geometry::Triangle g(*itF);
			area_minima += g.area();
		}
		area_minima /= 40000; // Questo limita a 10000 facce massime
	}
}

/*
 * Inserisce le condizioni al contorno
 */
void Poisson::setBConditions() {
	CDT::Vertex_circulator cV = cdt_.incident_vertices (cdt_.infinite_vertex());
	CDT::Vertex_circulator begin = cV;
	do {
		cV->info().toDirichlet();
		cV->info().value() = 0.;
		++cV;
	} while(cV != begin);
}

/*
 * Associa ad ogni nodo un indice
 */
void Poisson::enumNodes() {
	n_point = 0;
	CDT::Finite_vertices_iterator itV;
	for (itV = cdt_.finite_vertices_begin(); itV != cdt_.finite_vertices_end(); ++itV)
		if (!itV->info().isDirichlet())
			itV->info().index() = n_point++;
}

/*
 * Assemblaggio della matrice di stiffness e del termine noto
 */
void Poisson::makeMatrixTermineNoto(Matrix &A , Vector &F) {
	// Pattern della matrice sparsa
	{
		// Costruzione del pattern
		std::vector<std::list<unsigned int> > pattern;
		pattern.resize(n_point);
		unsigned int nnz = 0;
		{
			CDT::Finite_vertices_iterator itV;
			for (itV = cdt_.finite_vertices_begin(); itV != cdt_.finite_vertices_end(); ++itV) {
				if (!itV->info().isDirichlet()) {
					(pattern[itV->info().index()]).push_back(itV->info().index());
					++nnz;
					CDT::Vertex_circulator ciV = cdt_.incident_vertices(itV);
					CDT::Vertex_circulator end = ciV;
					do {
						if (!cdt_.is_infinite(ciV) && !ciV->info().isDirichlet()) {
							(pattern[itV->info().index()]).push_back(ciV->info().index());
							++nnz;
						}
						++ciV;
					} while (ciV != end);
				}
			}
		}
		// Inserimento del pattern
		A.newsize(n_point,n_point,nnz);
		{
			unsigned int i = 1;
			unsigned int j = 0;
			A.row_ptr(0) = 0;
			std::vector<std::list<unsigned int> >::const_iterator itv;
			for (itv = pattern.begin(); itv != pattern.end(); ++itv) {
				std::list<unsigned int>::const_iterator itl;
				A.row_ptr(i) = A.row_ptr(i-1) + itv->size();
				++i;
				for (itl = itv->begin(); itl != itv->end(); ++itl) {
					A.val(j) = 0.;
					A.col_ind(j) = *itl;
					++j;
				}
			}
		}
	}
	{
		F.newsize(n_point);
		F = 0.;
		CDT::Finite_faces_iterator itF;
		for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
			// Geometria dell'integratore e elemento
			Triangle g(*itF);
			ElementP1 e(g);
			Integrator2 integr(g);

			LinearAlgebra::Vector<unsigned int, ElementP1::DegreeOfFreedom> index;
			LinearAlgebra::Vector<bool, ElementP1::DegreeOfFreedom> cond;

			for (unsigned int i = 0; i < ElementP1::DegreeOfFreedom; ++i) {
				// Indici vertici del triangolo
				index(i) = itF->vertex(i)->info().index();
				// Condizioni al bordo
				cond(i) = !itF->vertex(i)->info().isDirichlet();
			}

			for (unsigned int i = 0; i < ElementP1::DegreeOfFreedom; ++i) {
				// Costruzione termine noto
				if (cond(i))
					F[index(i)] += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(forzante, e.base(i));
				// Termine di diffusione
				for (unsigned int j = 0; j < ElementP1::DegreeOfFreedom; ++j)
					if (cond(i) && cond(j))
						A.set(index(i),index(j)) += integr.integrateMul<Integrator2::Transform, Integrator2::Transform>(e.grad(i), e.grad(j));
			}
		}
	}
}

/*
 * Soluzione del sistema
 */
void Poisson::solveSystem(Matrix const & A, Vector & x, Vector const & b) {
	double tol = 1.e-6;
	int maxit = n_point;
	int result;

	x.newsize(n_point);
	x = 0.;
	
	DiagPreconditioner_double D(A);
	
	result = CG(A, x, b, D, maxit, tol);
}

/*
 * Salvataggio della soluzione
 */
void Poisson::saveSolution(Vector const & x) {
	CDT::Finite_vertices_iterator itV;
	for(itV = cdt_.finite_vertices_begin(); itV != cdt_.finite_vertices_end(); ++itV)
		if (!itV->info().isDirichlet())
			itV->info().value() = x[itV->info().index()];
}

/*
 * Raffinamento della mesh
 */
template < typename Stimator >
unsigned int Poisson::refine(double const & eps) {
	
	Stimator s(cdt_);
	
	double max_res = (3. * eps ) / (2. * std::sqrt(cdt_.number_of_faces ()));
	double min_res = max_res / 3.;
	
	PointList da_inserire;
	std::list<CDT::Vertex_handle> da_eliminare;
	
	CDT::Finite_faces_iterator itF;
	CDT::Finite_vertices_iterator itV;
	
	
	// Inserimento
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Calcolo residuo
		itF->info().res() = s.residue(itF, forzante);
		
		Geometry::Triangle g(*itF);
		if ((itF->info().res() > max_res) && (g.area() > area_minima)) {
			LinearAlgebra::Vector<double, 2> p0(itF->vertex(0)->point().x(), itF->vertex(0)->point().y());
			LinearAlgebra::Vector<double, 2> p1(itF->vertex(1)->point().x(), itF->vertex(1)->point().y());
			LinearAlgebra::Vector<double, 2> p2(itF->vertex(2)->point().x(), itF->vertex(2)->point().y());
			
			LinearAlgebra::Vector<double, 2> m((p0+p1+p2) / 3.);
			
			da_inserire.push_back(Point(m(0), m(1)));
		}
	}
	
	// Rimozione
	for (itV =cdt_.finite_vertices_begin(); itV != cdt_.finite_vertices_end(); ++itV) {
		// Non eliminabili
		if (!itV->info().isBoundary()) {
			// Residuo medio
			double t = 0.;
			unsigned int n = 0;
			CDT::Face_circulator cF = cdt_.incident_faces(itV);
			CDT::Face_circulator end = cF;
			do {
				if (!cdt_.is_infinite(cF)) {
					t += cF->info().res();
					++n;
				}
				++cF;
			} while(cF != end);
			t /= n;

			// Inserimento
			if (t < min_res) da_eliminare.push_back(itV);
		}
	}
	
	if (da_eliminare.size()) {
		std::list<CDT::Vertex_handle>::iterator it;
		for (it = da_eliminare.begin(); it != da_eliminare.end(); ++it)
			cdt_.remove(*it);
	}

	if (da_inserire.size())
		cdt_.insert(da_inserire.begin(), da_inserire.end());

	// Costruzione della griglia
	setBConditions();
	enumNodes();

	// Costruzione matrice termine noto
	Matrix A;
	Vector x;
	Vector F;	
	makeMatrixTermineNoto(A, F);

	// Soluzione del sistema lineare
	solveSystem(A, x, F);

	// Salvataggio della soluzione nella griglia
	saveSolution(x);

	return da_inserire.size() + da_eliminare.size();
}

/*
 * Calcolo dell'errore in norma L2
 */
template < typename Function >
double Poisson::errL2(Function const & f) {
	double err = 0.;
	CDT::Finite_faces_iterator itF;
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Soluzione
		LinearAlgebra::Vector<double, 3> u;
		for (unsigned int i = 0; i < 3; ++i)
			u(i) = itF->vertex(i)->info().value();

		// Elemento
		Geometry::Triangle g(*itF);
		ElementP1 e(g, u);

		double t = e.normL2(f);
		err += (t*t);
	}
	return std::sqrt(err);
}

/*
 * Calcolo dell'errore in norma H1
 */
template <typename Function1, typename Function2, typename Function3>
double Poisson::errH1(Function1 const & f, Function2 const & fx, Function3 const & fy) {
	double err = 0.;
	CDT::Finite_faces_iterator itF;
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Soluzione
		LinearAlgebra::Vector<double, 3> u;
		for (unsigned int i = 0; i < 3; ++i)
			u(i) = itF->vertex(i)->info().value();

		// Elemento
		Geometry::Triangle g(*itF);
		ElementP1 e(g, u);

		double t = e.normH1(f, fx, fy);
		err += (t*t);
	}
	return std::sqrt(err);
}

/*
 * Calcolo della norma L2
 */
double Poisson::normL2() {
	double n = 0.;
	CDT::Finite_faces_iterator itF;
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Soluzione
		LinearAlgebra::Vector<double, 3> u;
		for (unsigned int i = 0; i < 3; ++i)
			u(i) = itF->vertex(i)->info().value();

		// Elemento
		Geometry::Triangle g(*itF);
		ElementP1 e(g, u);

		double t = e.normL2();
		n += (t*t);
	}
	return std::sqrt(n);
}

/*
 * Calcolo della  norma H1
 */
double Poisson::normH1() {
	double n = 0.;
	CDT::Finite_faces_iterator itF;
	for (itF = cdt_.finite_faces_begin(); itF != cdt_.finite_faces_end(); ++itF) {
		// Soluzione
		LinearAlgebra::Vector<double, 3> u;
		for (unsigned int i = 0; i < 3; ++i)
			u(i) = itF->vertex(i)->info().value();

		// Elemento
		Geometry::Triangle g(*itF);
		ElementP1 e(g, u);

		double t = e.normH1();
		n += (t*t);
	}
	return std::sqrt(n);
}

#endif // POISSON_HPP
