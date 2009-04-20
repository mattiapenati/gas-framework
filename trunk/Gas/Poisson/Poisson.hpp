#ifndef POISSON_HPP
#define POISSON_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "Info.hpp"

#include "../SparseLib++/coord_double.h"
#include "../SparseLib++/comprow_double.h"
#include "../SparseLib++/diagpre_double.h"
#include "../Iml/cg.h"

#include "../Integration/Integration.h"

#include <list>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

class Poisson {
	private:
		// Kernel
		struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
		
		// Triangolazione
		typedef CGAL::Triangulation_vertex_base_with_info_2<PointInfo, K> Vb;
		typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
		
		typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
		
		typedef double (*func)(double const &, double const &);
	public:
		typedef CDT::Point Point;
		typedef std::list<Point> PointList;

	public:
		Poisson(PointList const &, func, double const & criteria = 0.);
		~Poisson();
		void saveToSVG(char const *, bool const);
	private:
		// LinearAlgebra
		typedef Coord_Mat_double Matrix;
		typedef MV_Vector_double Vector;
		
		// Integrator
		typedef Integrator < Method::NewtonCotes_2<Geometry::Triangle<CDT::Face>, 2> >  Integrator2;
		
		CDT cdt_;
		Integrator2 NC;
		
		unsigned int n_point;
		
		func dato_al_bordo;
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
		
};

/*
 * Costruttore di default
 */
Poisson::Poisson(PointList const &boundary, func f, double const & criteria) : NC() {
	
	dato_al_bordo = f;
	
	// Costruzione della griglia
	insertNodes(boundary, criteria);
	setBConditions();
	enumNodes();
	
	// Costruzione matrice termine noto
	Matrix A;
	Vector F;	
	makeMatrixTermineNoto(A, F);
	
	// Soluzione del sistema lineare
	Vector x(n_point, 0.);
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
	double xmin, xmax, ymin, ymax, umin, umax;
	double x, y, u;
	int color, r, g ,b;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	xmin = xmax = itV->point().x();
	ymin = ymax = itV->point().y();
	umin = umax = itV->info().value();
	++itV;
	while(itV != cdt_.finite_vertices_end()) {
		x = itV->point().x();
		y = itV->point().y();
		u = itV->info().value();
		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
		if (y < ymin) ymin = y;
		if (y > ymax) ymax = y;
		if (u > umax) umax = u;
		if (u < umin) umin = u;
		++itV;
	}
	W = H * (xmax -xmin) / (ymax - ymin);
	// Inizio SVG
	f << "<svg width=\""<<W<<"px\" height=\""<<H<<"px\" " \
	  << "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	// Triangolazione
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	while(itF != cdt_.finite_faces_end()) {
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
		
		u = 0.;
		u += itF->vertex(0)->info().value();
		u += itF->vertex(1)->info().value();
		u += itF->vertex(2)->info().value();
		u /= 3;
		
		// Calcolo dei colori
		color = 255 * (u - umin ) / (umax - umin);
		r = color;
		g = color;
		b = 0;
		
		// Stampa a file
		if (triangulation)
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"white\" width=\"1\">\n";
		else
			f << "<g id=\"triangulation\" fill=\"rgb("<<r<<", "<<g<<", "<<b<<")\" stroke=\"rgb("<<r<<", "<<g<<", "<<b<<")\" width=\"1\">\n";
		f << "\t<polygon points=\""<<static_cast<int>(x1)<<", "<<static_cast<int>(y1)<<" "\
			<< static_cast<int>(x2)<<", "<<static_cast<int>(y2)<<" "\
			<< static_cast<int>(x3)<<", "<<static_cast<int>(y3)<<"\" />\n";
		f << "</g>\n";
		++itF;
	}
	f  << "</svg>\n";
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
	
	x0 = (*itV1)->point().x();
	y0 = (*itV1)->point().y();
	x1 = (*itV2)->point().x();
	y1 = (*itV2)->point().y();
	
	minimo = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	
	while(itV2 != Pb_.end()) {
		cdt_.insert_constraint(*itV1, *itV2);
		
		// lunghezza del lato
		x0 = (*itV1)->point().x();
		y0 = (*itV1)->point().y();
		x1 = (*itV2)->point().x();
		y1 = (*itV2)->point().y();
	
		d = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
		if (d < minimo)
			minimo = d;
		
		// avanzamento
		++itV1;
		++itV2;
	}
	
	// Ultimo lato
	cdt_.insert_constraint(Pb_.front(), Pb_.back());
	
	x0 = Pb_.front()->point().x();
	y0 = Pb_.front()->point().y();
	x1 = Pb_.back()->point().x();
	y1 = Pb_.back()->point().y();

	d = std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	if (d < minimo)
		minimo = d;
	
	if (criteria != 0.)
		minimo = criteria;
	
	// Raffinamento
	CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, minimo));
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
	unsigned int i = 0u;
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	while(itV != cdt_.finite_vertices_end()) {
		if (!itV->info().isDirichlet()) {
			itV->info().index() = i;
			++i;
		}
		++itV;
	}
	n_point = i;
}

/*
 * Assemblaggio della matrice di stiffness e del termine noto
 */
void Poisson::makeMatrixTermineNoto(Matrix &A , Vector &F) {
	// Inserimento valori
	Vector Tmp;
	F.newsize(n_point);
	F = 0.;
	// TODO migliorare l'utilizzo delle matrici sparse
	Tmp.newsize(n_point * n_point);
	Tmp = 0.;
	CDT::Finite_faces_iterator itF = cdt_.finite_faces_begin();
	while(itF != cdt_.finite_faces_end()) {
		
		// Coordinate vertici del triangolo
		double x0 = itF->vertex(0)->point().x(); double y0 = itF->vertex(0)->point().y();
		double x1 = itF->vertex(1)->point().x(); double y1 = itF->vertex(1)->point().y();
		double x2 = itF->vertex(2)->point().x(); double y2 = itF->vertex(2)->point().y();
		
		// Geometria dell'integratore
		Integrator2::Geometry g(*itF);
		
		// Area del triangolo
		double detJ = g.det();
		
		// Questo serve per evitare un bug nella libreria CGAL
		detJ = std::abs(detJ);
		if (detJ <= 1.0e-4) {
			++itF;
			continue;
		}
		
		// Determinante dello Jacobiano
		double invdetJ = 1./detJ;
		
		// Indici vertici del triangolo
		int i0 = itF->vertex(0)->info().index();
		int i1 = itF->vertex(1)->info().index();
		int i2 = itF->vertex(2)->info().index();
		
		// Vettori ausiliari
		double phi01 = y0 - y2 + x2 - x0;
		double phi02 = y1 - y0 + x0 - x1;
		double phi11 = y2 - y0;
		double phi12 = y0 - y1;
		double phi21 = x0 - x2;
		double phi22 = x1 - x0;
		
		// Condizioni al bordo
		bool b0 = !itF->vertex(0)->info().isDirichlet();
		bool b1 = !itF->vertex(1)->info().isDirichlet();
		bool b2 = !itF->vertex(2)->info().isDirichlet();
		
		// Assemblaggio matrice di stiffness
		
		// Termine di diffusione
		if (b0 && b0) Tmp(n_point*i0+i0) += 0.5*invdetJ*((phi01*phi01) + (phi02*phi02));
		if (b0 && b1) Tmp(n_point*i0+i1) += 0.5*invdetJ*((phi11*phi01) + (phi12*phi02));
		if (b0 && b2) Tmp(n_point*i0+i2) += 0.5*invdetJ*((phi21*phi01) + (phi22*phi02));
		if (b1 && b0) Tmp(n_point*i1+i0) += 0.5*invdetJ*((phi01*phi11) + (phi02*phi12));
		if (b1 && b1) Tmp(n_point*i1+i1) += 0.5*invdetJ*((phi11*phi11) + (phi12*phi12));
		if (b1 && b2) Tmp(n_point*i1+i2) += 0.5*invdetJ*((phi21*phi11) + (phi22*phi12));
		if (b2 && b0) Tmp(n_point*i2+i0) += 0.5*invdetJ*((phi01*phi21) + (phi02*phi22));
		if (b2 && b1) Tmp(n_point*i2+i1) += 0.5*invdetJ*((phi11*phi21) + (phi12*phi22));
		if (b2 && b2) Tmp(n_point*i2+i2) += 0.5*invdetJ*((phi21*phi21) + (phi22*phi22));
		
		// Definizione del dominio dell'integrazione
		NC.domain(g);
		
		// Costruzione termine noto
		if (b0) F[i0] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(dato_al_bordo, Poisson::P1_phi0);
		if (b1) F[i1] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(dato_al_bordo, Poisson::P1_phi1);
		if (b2) F[i2] += NC.integrateMul<Integrator2::Transform, Integrator2::NoTransform>(dato_al_bordo, Poisson::P1_phi2);
		
		// Avanzamento dell'iteratore
		++itF;
	}
	unsigned int nz = 0;
	for(int i = 0; i < n_point*n_point; ++i) {
		if (Tmp(i) != 0.)
			++nz;
	}
	A.newsize(n_point, n_point, nz);
	// Inserimento nella matrice
	int j = 0;
	for(int i = 0; i < n_point*n_point; ++i) {
		if(Tmp(i) != 0.) {
			// Valore
			A.val(j) = Tmp(i);
			// Indice colonna
			A.col_ind(j) = i % n_point;
			// Indice riga
			A.row_ind(j) = i / n_point;
			j++;
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
	CompRow_Mat_double Atmp(A);
	
	DiagPreconditioner_double D(Atmp);
	
	result = CG(Atmp, x, b, D, maxit, tol);
}

/*
 * Salvataggio della soluzione
 */
void Poisson::saveSolution(Vector const & x) {
	CDT::Finite_vertices_iterator itV = cdt_.finite_vertices_begin();
	while(itV != cdt_.finite_vertices_end()) {
		if (!itV->info().isDirichlet()) {
			itV->info().value() = x[itV->info().index()];
		}
		++itV;
	}
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
				pbm<<"1 ";
			else
				pbm<<"0 ";
			pbm<<std::endl;
		}
	}
}

#endif // POISSON_HPP